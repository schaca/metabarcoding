#!/bin/bash

echo "fastqc"
~/Programs/FastQC/fastqc COI/0_raw/*.fastq -o COI/1_fastqc/

echo "merge"
for file in ./COI/0_raw/*_R1_001.fastq ; do
    filename=$(basename "$file")
    illuminapairedend --score-min=40 -r ./COI/0_raw/$(echo "$filename" | cut -d "_" -f1,2,3)_R2_001.fastq $file > COI/2_merge/$(echo "$filename" | cut -d "_" -f1,2).merged.fastq
done

for file in ./COI/2_merge/*.fastq ; do
    filename=$(basename "$file" | cut -d "." -f1)
    echo "obigrep"
    obigrep -p 'mode!="joined"' $file > ./COI/3_align_filter/$filename.merged.ali.fastq
    echo "tagcleaner"
    perl ~/Programs/tagcleaner-standalone-0.16/tagcleaner.pl -fastq ./COI/3_align_filter/$filename.merged.ali.fastq -out ./COI/4_primerremoval/$filename.merged.ali.pr -tag5 GGWACWGGWTGAACWGTWTAYCCYCC -tag3 TGRTTYTTYGGCCAYCCCGARGTCTA -mm5 3 -mm3 3 -info
    echo "obiuniq"
    obiuniq -m sample ./COI/4_primerremoval/$filename.merged.ali.pr.fastq > ./COI/5_unique+ann_by_count/$filename.merged.ali.pr.uniq.fasta  
    echo "obiannotate"    
    obiannotate -k count ./COI/5_unique+ann_by_count/$filename.merged.ali.pr.uniq.fasta > $$ ; mv $$ ./COI/5_unique+ann_by_count/$filename.merged.ali.pr.uniq.fasta
    echo "obigrep"    
    obigrep -l 200 -L 400 -p 'count>=10' ./COI/5_unique+ann_by_count/$filename.merged.ali.pr.uniq.fasta > ./COI/6_filter/$filename.merged.ali.pr.uniq.c10.l200.L400.fasta
    echo "obiclean"
    obiclean -r 0.05 -H ./COI/6_filter/$filename.merged.ali.pr.uniq.c10.l200.L400.fasta > ./COI/7_clean/$filename.merged.ali.pr.uniq.c10.l200.L400.clean.fasta
    echo "blast"
    ~/Programs/ncbi-blast-2.5.0+/bin/blastn -query ./COI/7_clean/$filename.merged.ali.pr.uniq.c10.l200.L400.clean.fasta -db nt -remote -max_target_seqs 1 -outfmt "6 qseqid stitle length pident" -out ./COI/8_blastn/$filename.identifications
done
