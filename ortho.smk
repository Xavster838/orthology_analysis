import os
SMK_DIR = os.path.dirname(workflow.snakefile) #workflow.snakefile == full path to snakefile 
shell.prefix("source {}/env.cfg; set -eo pipefail".format(SMK_DIR)) #source env config before anything gets run ##these eo sets tells pipeline to stop running if one of lines fails.
configfile: "ortho.yaml"
workdir: "ortho_results" #add output directory
SMS = list(config.keys()) #just strings of element names in yaml ["a","b"]
SM_region = {}
SM_qs = {}
SM_ref = {}
for SM in SMS:
    SM_region[SM] = config[SM]["regions"]
    SM_qs[SM] = config[SM]["qs"]
    SM_ref[SM] = config[SM]["ref"]

wildcard_constraints:
    SM = "|".join(SMS)


def get_ref(wc): #wc is snakemake object that looks at wildcards of rule that calls get_ref
    return( SM_ref[wc.SM] )
def get_qs(wc):
    return( SM_qs[wc.SM])

rule all:
    input:
        bams = expand("{SM}.bam", SM = SMS) 
        

rule alignment:
    input:
        ref = get_ref , 
        qs = get_qs ,
    output:
        bam = "{SM}.bam"
    threads: 16 #snakemake will scale to either 16 or largest n_threads it can.
    shell:""" 
echo {input}
minimap2 -x asm20 -a --eqx -Y -r 50000 -M 0 --hard-mask-level --secondary=no \
        -t {threads} -I 8g -2 -K 1500m \
        {input.ref} <(cat {input.qs}) \
        | samtools view -b - | samtools sort - > {output.bam}
    """

