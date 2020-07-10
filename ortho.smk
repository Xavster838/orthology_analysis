import os
SMK_DIR = os.path.dirname(workflow.snakefile) #workflow.snakefile == full path to snakefile 
shell.prefix("source {}/env.cfg; set -eo pipefail;".format(SMK_DIR)) #source env config before anything gets run ##these eo sets tells pipeline to stop running if one of lines fails.
configfile: "ortho.yaml"
SMS = list(config.keys()) #just strings of element names in yaml ["a","b"]
SM_region = {}
SM_qs = {}
SM_ref = {}

def get_regions_list( bed_file_path ):
    """
    given bed_file_path return list of region element characters.
    """
    def get_region(line):
        """
        helper function: given line from bed file, process and turn into regions format.
        """
        line_split = line.strip().split() #strip flanking whitespace and split on whitespace
        return("{}:{}-{}".format(line_split[0], line_split[1], line_split[2] ))
    
    regions_list = []
    with open(bed_file_path) as fp:
        for cnt, line in enumerate(fp):
            if(line[0] == '#'): #skip comments in bed file.
                continue
            regions_list.append( get_region( str(line) ) )
    
    return(regions_list)


for SM in SMS:
    SM_region[SM] = get_regions_list(config[SM]["regions"])
    SM_qs[SM] = [os.path.abspath(q) for q in config[SM]["qs"] ]
    SM_ref[SM] = os.path.abspath(config[SM]["ref"])

workdir: "ortho_results" #add output directory

wildcard_constraints:
    SM = "|".join(SMS)


def get_ref(wc): #wc is snakemake object that looks at wildcards of rule that calls get_ref
    return( SM_ref[wc.SM] )
def get_qs(wc):
    return( SM_qs[wc.SM])

def get_region_fastas(wc):
    all_fastas = []
    for SM in SMS:
        cur_regions = SM_region[SM]

        fastas = expand("{SM}_{RGN}.fa", SM = [SM] , RGN = range( len(cur_regions) ) ) #double expand...
        all_fastas += fastas
    return( all_fastas )



rule all:
    input:
        fastas = get_region_fastas 
    
rule alignment:
    input:
        ref = get_ref , 
        qs = get_qs ,
    output:
        bam = "{SM}.bam",
    threads: 16 #snakemake will scale to either 16 or largest n_threads it can.
    shell:""" 
echo {input}
minimap2 -x asm20 -a --eqx -Y -r 50000 -M 0 --hard-mask-level --secondary=no \
        -t {threads} -I 8g -2 -K 1500m \
        {input.ref} <(cat {input.qs}) \
        | samtools view -b - | samtools sort - > {output.bam}
    """

rule get_alignment_index:
    input:
        bam = rules.alignment.output.bam
    output:
        bai = rules.alignment.output.bam + ".bai"
    shell:"""
module list
samtools index {input.bam}
"""

def get_region( wc  ):
    return SM_region[wc.SM][int(wc.RGN)]

rule region_fasta:
    input:
        bam = "{SM}.bam",
        bai = rules.get_alignment_index.output.bai ,
        fasta = get_ref 
    params:
        regions = get_region 
    output:
        fasta = "{SM}_{RGN}.fa" #note may have to change the naming scheme...
    shell:"""
#echo "{params.regions}"
subseq_path=/net/eichler/vol27/projects/structural_variation/nobackups/tools/seqtools/201910/CentOS6/bin/subseqfa
$subseq_path -r '{params.regions}' -b -o {output.fasta} {input.bam}

samtools faidx {input.fasta} '{params.regions}' >> {output.fasta}
"""






