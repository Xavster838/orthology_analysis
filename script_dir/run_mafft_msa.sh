#!/bin/bash

INPUT_FASTA=$1
OUTPUT_FASTA_NAME=$2

#test that there is more than 1 sequence
n_seq=`grep ">" -c ${INPUT_FASTA}`

if [[ $n_seq -gt 1 ]]
then
    mafft --reorder --treeout $INPUT_FASTA > ${OUTPUT_FASTA_NAME}
else
    touch ${OUTPUT_FASTA_NAME}
fi
