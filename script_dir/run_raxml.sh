#!/bin/bash
# run raxml in MSA step for given regions.
# inputs: output directory, MSA file .
# output: all outputs of raxml pipeline.
INPUT_MSA_fasta=$1
OUTPUT_name=$2
OUTPUT_dir=$3


initial_dir=$PWD
script_dir=`dirname $0`
input_dir=`dirname ${INPUT_MSA_fasta}`


echo $script_dir
echo `ls ${script_dir}`
echo $input_dir
echo `ls ${input_dir}`

echo "switching to ${OUTPUT_dir}"
mkdir -p ${OUTPUT_dir}
cd ${OUTPUT_dir}

if [ -s ${INPUT_MSA_fasta} ]
then
    echo "running raxml."
    echo "running on ${INPUT_MSA_fasta}."
    echo "editing sequence headers to remove : character."
    sed -e 's/:/–/g' ${INPUT_MSA_fasta} > temp_corrected_name #${INPUT_MSA_fasta}_temp_character_corrected.fa 
    echo "converting to phylum format."
    ${script_dir}/fasta_to_phy.sh temp_corrected_name > temp_corrected_name.phy
    echo "running raxml."
    ~pdx/bin/raxml -p 12345 -x 12345 -s temp_corrected_name.phy -m GTRGAMMA -# 100 -T 8 -n ${OUTPUT_name} 
    #make sure output file made.
    touch RAxML_bootstrap.${OUTPUT_name}
else
    touch RAxML_bootstrap.${OUTPUT_name}  
fi

echo "Done. Returning to original directory."
cd $initial_dir
