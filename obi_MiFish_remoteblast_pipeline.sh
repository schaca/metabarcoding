#!/bin/bash

IN_DIR="./raw/"
OUT_DIR="./obi_remoteblast"
rm -rf $OUT_DIR
mkdir $OUT_DIR
mkdir $OUT_DIR/0_fastqc
mkdir $OUT_DIR/1_merge
mkdir $OUT_DIR/2_filter
mkdir $OUT_DIR/3_primerremoval
mkdir $OUT_DIR/4_unique+ann-by-count
mkdir $OUT_DIR/5_filter
mkdir $OUT_DIR/6_clean
mkdir $OUT_DIR/7_blastn
mkdir $OUT_DIR/8_results

echo -e "Output folders created > $OUT_DIR/"

#--------------#
for i in $(seq 1 $(stty size | cut -d ' ' -f2)); do 
    echo -n "-" 
done
#--------------#

for file in $IN_DIR/*R1_001.fastq ; do
    filename=$(basename "$file" | cut -d "_" -f1,2,3)
    samplename=$(basename "$file" | cut -d "-" -f1)

    echo -e "Processing starts for sample : $samplename"

    echo -e "FastQC(v0.11.5) starts"
    ~/Programs/FastQC/fastqc $IN_DIR/"$filename"_R1_001.fastq -o $OUT_DIR/0_fastqc/
    ~/Programs/FastQC/fastqc $IN_DIR/"$filename"_R2_001.fastq -o $OUT_DIR/0_fastqc/
    echo -e "FastQC results recorded > $OUT_DIR/0_fastqc/"

    echo -e "Merging by OBITOOLS-illuminapairedend(min.ali.score=40)"
    illuminapairedend --score-min=40 -r $IN_DIR/"$filename"_R2_001.fastq $IN_DIR/"$filename"_R1_001.fastq > $OUT_DIR/1_merge/$samplename.merged.fastq
    echo -e "Merged sequences are recorded > $OUT_DIR/1_merge/"

    echo "Filtering unmerged sequences by OBITOOLS-obigrep"
    obigrep -p 'mode!="joined"' $OUT_DIR/1_merge/$samplename.merged.fastq > $OUT_DIR/2_filter/$samplename.aligned.fastq
    echo -e "Filtered sequences are recorded > $OUT_DIR/2_filter/"

    echo "Trimming primer sequences by tagcleaner(v0.16)(3 mismatches allowed)"
    perl ~/Programs/tagcleaner-standalone-0.16/tagcleaner.pl -fastq $OUT_DIR/2_filter/$samplename.aligned.fastq -out $OUT_DIR/3_primerremoval/$samplename.aligned.trimmed -tag5 GTCGGTAAAACTCGTGCCAGC -tag3 CAAACTGGGATTAGATACCCCACTATG -mm5 3 -mm3 3 -info
    echo -e "Trimmed sequences are recorded > $OUT_DIR/3_primerremoval/"

    echo "Dereplicating sequences by OBITOOLS-obiuniq"
    obiuniq -m sample $OUT_DIR/3_primerremoval/$samplename.aligned.trimmed.fastq > $OUT_DIR/4_unique+ann-by-count/$samplename.aligned.trimmed.uniq.fasta
    echo "Annotation by count by OBITOOLS-obiannotate"
    obiannotate -k count $OUT_DIR/4_unique+ann-by-count/$samplename.aligned.trimmed.uniq.fasta > $$ ; mv $$ $OUT_DIR/4_unique+ann-by-count/$samplename.aligned.trimmed.uniq.fasta
    echo -e "Dereplicated and annotated sequences are recorded > $OUT_DIR/4_unique+ann-by-count"    

    echo "Length and count filtering by OBITOOLS-obigrep(minlen=100,maxlen=300)"
    echo "Filtering counts<10" 
    obigrep -l 100 -L 300 -p 'count>=10' $OUT_DIR/4_unique+ann-by-count/$samplename.aligned.trimmed.uniq.fasta > $OUT_DIR/5_filter/$samplename.aligned.trimmed.uniq.c10.l100.L300.fasta
    echo "Filtering counts >10 and <5" 
    obigrep -l 100 -L 300 -p 'count<10' $OUT_DIR/4_unique+ann-by-count/$samplename.aligned.trimmed.uniq.fasta | obigrep -p 'count>=5' > $OUT_DIR/5_filter/$samplename.aligned.trimmed.uniq.c5-10.l100.L300.fasta
    echo "Filtering counts >5"
    obigrep -l 100 -L 300 -p 'count<5' $OUT_DIR/4_unique+ann-by-count/$samplename.aligned.trimmed.uniq.fasta | obigrep -p 'count>=3' > $OUT_DIR/5_filter/$samplename.aligned.trimmed.uniq.c3-5.l100.L300.fasta
    echo -e "Quality filtered sequences are recorded > $OUT_DIR/5_filter"    
    
    echo "Filtering PCR/sequencing errors by OBITOOLS-obiclean(r=0.05)"
    obiclean -r 0.05 -H $OUT_DIR/5_filter/$samplename.aligned.trimmed.uniq.c10.l100.L300.fasta > $OUT_DIR/6_clean/$samplename.aligned.trimmed.uniq.c10.l100.L300.clean.fasta
    obiclean -r 0.05 -H $OUT_DIR/5_filter/$samplename.aligned.trimmed.uniq.c5-10.l100.L300.fasta > $OUT_DIR/6_clean/$samplename.aligned.trimmed.uniq.c5-10.l100.L300.clean.fasta
    obiclean -r 0.05 -H $OUT_DIR/5_filter/$samplename.aligned.trimmed.uniq.c3-5.l100.L300.fasta > $OUT_DIR/6_clean/$samplename.aligned.trimmed.uniq.c3-5.l100.L300.clean.fasta
    echo -e "Filtered sequences are recorded > $OUT_DIR/6_clean"
    
    echo "Local BLAST (ncbi-blast-2.5.0+) starts"
    ~/Programs/ncbi-blast-2.5.0+/bin/blastn -query $OUT_DIR/6_clean/$samplename.aligned.trimmed.uniq.c10.l100.L300.clean.fasta -db nt -remote -max_target_seqs 1 -outfmt "6 qseqid stitle length pident" -out $OUT_DIR/7_blastn/$samplename.c10.identifications
    ~/Programs/ncbi-blast-2.5.0+/bin/blastn -query $OUT_DIR/6_clean/$samplename.aligned.trimmed.uniq.c5-10.l100.L300.clean.fasta -db nt -remote -max_target_seqs 1 -outfmt "6 qseqid stitle length pident" -out $OUT_DIR/7_blastn/$samplename.c5-10.identifications
    ~/Programs/ncbi-blast-2.5.0+/bin/blastn -query $OUT_DIR/6_clean/$samplename.aligned.trimmed.uniq.c3-5.l100.L300.clean.fasta -db nt -remote -max_target_seqs 1 -outfmt "6 qseqid stitle length pident" -out $OUT_DIR/7_blastn/$samplename.c3-5.identifications
    echo -e "Identifications are recorded > $OUT_DIR/7_blastn"

    echo "Parsing BLAST results"
    bash ~/Programs/Scripts/obi_blast_parser.sh
    
    echo "All results recorded > $OUT_DIR/8_results"
done

echo -e "FINISHED!!"
