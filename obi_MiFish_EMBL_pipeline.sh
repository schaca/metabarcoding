#!/bin/bash

for file in ./*R1_001.fastq ; do
    filename=$(basename "$file" | cut -d "_" -f1,2,3)
    illuminapairedend --score-min=40 -r "$filename"_R2_001.fastq "$filename"_R1_001.fastq > "$filename".merged.fastq
done

for file in ./*_MiFish_L001.merged.fastq ; do
    filename=$(basename "$file" | cut -d "_" -f1,2)
    echo "obigrep"
    obigrep -p 'mode!="joined"' $file > $filename.ali.fastq
    echo "tagcleaner"
    perl ~/Programs/tagcleaner-standalone-0.16/tagcleaner.pl -fastq $filename.ali.fastq -out $filename.ali.pr -tag5 GTCGGTAAAACTCGTGCCAGC -tag3 CAAACTGGGATTAGATACCCCACTATG -mm5 3 -mm3 3 -info
    echo "obiuniq"
    obiuniq -m sample $filename.ali.pr.fastq > $filename.ali.pr.uniq.fasta 
    echo "obiann"    
    obiannotate -k count $filename.ali.pr.uniq.fasta > $$ ; mv $$ $filename.ali.pr.uniq.fasta 
    echo "obigrep"    
    obigrep -l 100 -L 300 -p 'count>=10' $filename.ali.pr.uniq.fasta > $filename.ali.pr.uniq.c10.l100.L300.fasta
    echo "obiclean"    
    obiclean -r 0.05 -H $filename.ali.pr.uniq.c10.l100.L300.fasta > $filename.ali.pr.uniq.c10.l100.L300.clean.fasta 
    echo "ecotag"    
    ecotag -d ~/Documents/Databases/EMBL/embl_vrt_r130 -R ~/Documents/Databases/EMBL/embl_vrt_r130.ref.fasta $filename.ali.pr.uniq.c10.l100.L300.clean.fasta > $filename.ali.pr.uniq.c10.l100.L300.clean.tag.fasta
    echo "obiann"    
    obiannotate  --delete-tag=scientific_name_by_db --delete-tag=obiclean_samplecount --delete-tag=obiclean_count --delete-tag=obiclean_singletoncount --delete-tag=obiclean_cluster --delete-tag=obiclean_internalcount --delete-tag=obiclean_head --delete-tag=taxid_by_db --delete-tag=obiclean_headcount --delete-tag=id_status --delete-tag=rank_by_db --delete-tag=order_name --delete-tag=order $filename.ali.pr.uniq.c10.l100.L300.clean.tag.fasta >  $filename.ali.pr.uniq.c10.l100.L300.clean.tag.ann.fasta
    echo "obisort"    
    obisort -k count -r $filename.ali.pr.uniq.c10.l100.L300.clean.tag.ann.fasta > $filename.ali.pr.uniq.c10.l100.L300.clean.tag.ann.sort.fasta
    echo "obisort"    
    obitab -o $filename.ali.pr.uniq.c10.l100.L300.clean.tag.ann.sort.fasta >  $filename.ali.pr.uniq.c10.l100.L300.clean.tag.ann.sort.tab
done

