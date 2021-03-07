#!/bin/bash

OUT_DIR="./"

for file in $OUT_DIR/7_blastn/*.identifications; do
    samplename=$(basename "$file" | cut -d "." -f1) 
    countnumber=$(basename "$file" | cut -d "." -f2)
    [ -d $OUT_DIR/8_results/$samplename ] || mkdir $OUT_DIR/8_results/$samplename        
    IFS=$'\n'
    for line in $(cut -f1 $file); do
        grep $line $OUT_DIR/6_clean/$samplename.aligned.trimmed.uniq.$countnumber.l100.L300.clean.fasta | cut -d " " -f2 | awk -v FS="(=|;)" '{print $2}' >> $OUT_DIR/8_results/$samplename/counts
    done
    paste $file $OUT_DIR/8_results/$samplename/counts | sort -t$'\t' -nrk2 > $OUT_DIR/8_results/$samplename/$samplename.$countnumber.table
    rm $OUT_DIR/8_results/$samplename/counts
    awk -F "\t" '$4>=98' $OUT_DIR/8_results/$samplename/$samplename.$countnumber.table > $OUT_DIR/8_results/$samplename/$samplename.$countnumber.98-100.table
    awk -F "\t" '$4<98' $OUT_DIR/8_results/$samplename/$samplename.$countnumber.table |  awk -F "\t" '$4>=95' > $OUT_DIR/8_results/$samplename/$samplename.$countnumber.95-98.table
    cut -f2 $OUT_DIR/8_results/$samplename/$samplename.$countnumber.98-100.table | cut -d " " -f1,2 | sort | uniq > $OUT_DIR/8_results/$samplename/$samplename.$countnumber.98-100.splist
    cut -f2 $OUT_DIR/8_results/$samplename/$samplename.$countnumber.95-98.table | cut -d " " -f1,2 | sort | uniq > $OUT_DIR/8_results/$samplename/$samplename.$countnumber.95-98.splist
        
    IFS=$'\n'
    for line in $(cat $OUT_DIR/8_results/$samplename/$samplename.$countnumber.98-100.splist); do
        grep $line $OUT_DIR/8_results/$samplename/$samplename.$countnumber.98-100.table | cut -f5 | awk '{ sum += $1 } END { print sum }' >> $OUT_DIR/8_results/$samplename/totalcounts
    done
    paste $OUT_DIR/8_results/$samplename/$samplename.$countnumber.98-100.splist $OUT_DIR/8_results/$samplename/totalcounts | sort -t$'\t' -nrk2 > $OUT_DIR/8_results/$samplename/$samplename.$countnumber.98-100.counts
    rm $OUT_DIR/8_results/$samplename/totalcounts

        
    IFS=$'\n'
    for line in $(cat $OUT_DIR/8_results/$samplename/$samplename.$countnumber.95-98.splist); do
        grep $line $OUT_DIR/8_results/$samplename/$samplename.$countnumber.95-98.table | cut -f5 | awk '{ sum += $1 } END { print sum }' >> $OUT_DIR/8_results/$samplename/totalcounts
    done
    paste $OUT_DIR/8_results/$samplename/$samplename.$countnumber.95-98.splist $OUT_DIR/8_results/$samplename/totalcounts | sort -t$'\t' -nrk2 > $OUT_DIR/8_results/$samplename/$samplename.$countnumber.95-98.counts
    rm $OUT_DIR/8_results/$samplename/totalcounts
    rm $OUT_DIR/8_results/$samplename/$samplename.$countnumber.98-100.splist
    rm $OUT_DIR/8_results/$samplename/$samplename.$countnumber.95-98.splist

done
