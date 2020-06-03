#!/bin/bash

files=""
if [[ $# == 1 && -d $1 ]]; then
    files=`find $1 -name "*.log"`
else
    files=$@
fi

grep -B 1 "sample time" $files |
    awk 'BEGIN {
           sampleType[0]="leaf"
           sampleType[1]="subtree"
           sampleType[2]="hap"
         }
         $0 ~ /resample_arg_regions/ {subtree=1};
         $0 ~ /resample_arg_leaf/ {subtree=0}
         $0 ~ /resample_arg_by_hap/ {subtree=2}
         {val=-1}
         $NF=="ms" {
            val=$(NF-1)/1000
         }
         $NF == "s" {
            val = $(NF-1)
         }
         $NF == "m" {
            val = $(NF-1)*60
         }
         $NF == "h" {
            val = $(NF-1)*3600
         }
         val > 0 {
            totalSec[subtree] += val;
            totalCount[subtree]++;
         }
         END{
            for (i=0; i < 3; i++) {
               if (totalCount[i] > 0) {
                  avg=totalSec[i]/totalCount[i]/60
                  print sampleType[i]": "avg" m (out of "totalCount[i]")"
                }
             }
         }'
