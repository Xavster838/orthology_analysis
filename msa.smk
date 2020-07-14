import os
SMK_DIR = os.path.dirname(workflow.snakefile) #workflow.snakefile == full path to snakefile 
SCRIPT_DIR = "{SMK_DIR}/script_dir"
shell.prefix("source {}/env.cfg; set -eo pipefail;".format(SMK_DIR)) #source env config before anything gets run ##these eo sets tells pipeline to stop running if one of lines fails.
configfile: "msa.yaml"
fofn = [os.path.abspath(line.strip()) for line in open(config["fofn"]) ] #get all fastas in fofn

def get_region_sample_name(region_fasta):
    """ get the name of the fasta file without the suffix for later running
    """
    base=os.path.basename(region_fasta)
    return os.path.splitext(base)[0]

region_2_fasta = {}
for fasta in fofn:
    key = get_region_sample_name(fasta)
    region_2_fasta[key] = fasta

RGNS = region_2_fasta.keys()

workdir: "MSA_results"

#wildcard_constraints:
#    fastas = "|".join([ get_region_sample_name(fasta) for fasta in config_fastas]) 
#RAxML_bootstrap.hg38_tb1d3_macaque_outgroup
def get_final_output_names(wc):
    """function to get final output names from initial input names"""
    output_names = []
    for fasta in config_fastas:
        cur_out = "{RGN}/RAxML_bootstrap.{RGN}".format(RGN = get_region_sample_name(fasta) )
        output_names.append(cur_out)
    return(output_names)

rule all:
    input:
        raxml = expand("{RGN}/RAxML_bootstrap.{RGN}", RGN = RGNS) 

def get_region_fasta(wc):
    return( region_2_fasta[wc.RGN] )

rule run_MSA:
    input:
        fasta = get_region_fasta 
    output:
        msa = "{RGN}/{RGN}_MSA.fa"
    shell:''''''
rule run_raxml:
    input:
        msa = rules.run_MSA.output.msa #"{RGN}/{RGN}_MSA.fa"
    output:
        raxml = "{RGN}/RAxML_bootstrap.{RGN}" 
    shell: '''
echo "working on {wildcards.RGN} snakemake" '''


       
