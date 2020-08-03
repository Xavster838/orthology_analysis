import os
SMK_DIR = os.path.dirname(workflow.snakefile) #workflow.snakefile == full path to snakefile 
shell.prefix("source {}/env.cfg; set -eo pipefail;".format(SMK_DIR)) #source env config before anything gets run ##these eo sets tells pipeline to stop running if one of lines fails.
configfile: "ortho.yaml"
CNF_DIR = os.path.dirname( os.path.abspath("ortho.yaml") )
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
		#line_split = line.strip().split() #strip flanking whitespace and split on whitespace
		return("{}:{}-{}".format(line[0], line[1], line[2] ))
    
	regions_dict = {} #regions_list = []
	
	for line in open(bed_file_path):
		if(line[0] == '#' or not line.strip()): #remove initial bed comment lines or empty lines
			continue
		t = line.strip().split() #strip white space and split into array of bed columns
		if( len(t) != 4):
			key = "_".join(t)
		else:
			key = t[3]
		regions_dict.setdefault(key , [])
		regions_dict[key] = get_region(t)
	#print(regions_dict.items())
	return(regions_dict)#return(regions_list)


for SM in SMS:
	SM_region[SM] = get_regions_list(config[SM]["regions"])
	SM_qs[SM] = [os.path.abspath(q) for q in config[SM]["qs"] ]
	SM_ref[SM] = os.path.abspath(config[SM]["ref"])


workdir: "ortho_results" #add output directory

wildcard_constraints:
	SM = "|".join(SMS) ,
	RGN = "|".join( [key for k , v in SM_region.items() for key in v]  ) #"\d+"
#ruleorder: region_fasta > query_region_fasta

def get_ref(wc): #wc is snakemake object that looks at wildcards of rule that calls get_ref
	return( SM_ref[wc.SM] )
def get_qs(wc):
	return( SM_qs[wc.SM])
def get_q_contigs( query_list):
	"""
	return list of all contigs in query list files. 
	"""
	contigs = []
	for cur_q in query_list:
		with open(cur_q) as fp:
			for cnt , line in enumerate(fp):
				if(line[0] == ">"):
					contigs.append( line.strip(">").rstrip().split()[0]) #get rif of > and ending white space, take first part of header. 
	return contigs

def get_region_fastas(wc):
    all_fastas = []
    for SM in SMS:
        cur_regions = SM_region[SM]

        fastas = expand("{SM}_{RGN}.fa", SM = [SM] , RGN = list(cur_regions.keys()) ) #double expand...
        all_fastas += fastas
    return( all_fastas )

def get_region_q_contig_sams(wc):
	"""
	get file names for fastas with region and query contig specific sequence.
	@input: wc
	@output: list of fasta file names used for rule all
	"""
	all_sams = []
	for SM in SMS:
		cur_regions = SM_region[SM]
		cur_q_contigs = get_q_contigs(SM_qs[SM])
		cur_SM_RGNs = expand("{SM}_{RGN}", SM = [SM] , RGN = list(cur_regions.keys() ) )  #double expand...
		for SM_RGN in cur_SM_RGNs:
			sams = expand(SM_RGN + "_{Q_CONTIG}.sam", Q_CONTIG = cur_q_contigs )
			all_sams += sams
	return( all_sams )

rule all:
	input:
		table = expand("{SM}_summary.tbl", SM = SMS) 
		#rgn_fastas = get_region_fastas ,
		#bed = get_perID_beds 

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
    return SM_region[wc.SM][wc.RGN]

rule region_fasta:
    input:
        bam = "{SM}.bam",
        bai = rules.get_alignment_index.output.bai ,
        fasta = get_ref 
    params:
        regions = get_region 
    output:
        fasta = "{SM}_{RGN}.fa" , 
        tmp = temp("tmp_{SM}_{RGN}.fa")
    shell:"""
#run subseq
subseq_path=/net/eichler/vol27/projects/structural_variation/nobackups/tools/seqtools/201910/CentOS6/bin/subseqfa
$subseq_path -r '{params.regions}' -b -o {output.tmp} {input.bam}

#get reference sequence and add to output fasta
samtools faidx {input.fasta} '{params.regions}' >> {output.tmp}

#fix name here with awk
echo {wildcards.RGN}
awk -v rgn={wildcards.RGN} '{{if ($0 ~ />/) {{printf "%s_%s",$1, rgn; print"" }} else {{print $0}} }}' {output.tmp} > {output.fasta}
"""

rule query_region_bam:
	input:
		ref = get_ref ,
		q_fasta = rules.region_fasta.output.fasta ,#"{SM}_{RGN}.fa" 
		#bam = "{SM}.bam" ,
		#bai = rules.get_alignment_index.output.bai
	params:
		regions = get_region
	output:
		query_region_bam = temp("query_region_bam/{SM}_{RGN}_{Q}.bam") ,
		bai = temp("query_region_bam/{SM}_{RGN}_{Q}.bam.bai")
	shell:"""
#run minimap2
minimap2 -x asm20 -a --eqx -Y -r 50000 -M 0 --hard-mask-level --secondary=no \
        -t {threads} -I 8g -2 -K 1500m \
        {input.ref} <(cat {input.q_fasta}) \
        | samtools view -b -S - | samtools view -b -F 2308 | samtools sort - > {output.query_region_bam}

#sort and index
samtools index {output.query_region_bam} 
"""

rule get_perID_results:
	input:
		bam_region_q = rules.query_region_bam.output.query_region_bam 
	output:
		bed = temp("perID/{SM}_{RGN}_{Q}.bed")
	shell:"""
perID_path=~mvollger/projects/hifi_asm/scripts/samIdentity.py
python3 $perID_path --bed {input.bam_region_q} > {output.bed}
"""

def get_perID_beds(wc):
    """
    get file names for beds processed by mvolleger's samIdentity.py.
    @input: wc
    @output: list of bed file names used for rule all
    """
    all_beds = []
    for SM in SMS:
        cur_regions = SM_region[SM]
        cur_q_contigs = get_q_contigs(SM_qs[SM])
        cur_SM_RGNs = expand("{SM}_{RGN}", SM = [SM] , RGN = list(cur_regions.keys() ) )  #double expand...
        for SM_RGN in cur_SM_RGNs:
            beds = expand(SM_RGN + "_{Q_CONTIG}.bed", Q_CONTIG = cur_q_contigs )
            all_beds += beds
    return( expand("perID/" + "{perID}" , perID = all_beds)) #all_beds )

def get_rgn_bed(wc):
	return "{dir_path}/{file_name}".format(dir_path = CNF_DIR, file_name = config[wc.SM]["regions"])

rule consolidate_table:
	input:
		perID = get_perID_beds ,
		rgn_bed = get_rgn_bed
	output:
		table = "{SM}_summary.tbl"
	shell:"""
in_dir=`dirname {input.perID[0]}`
echo ${{in_dir}}
files=(`ls -d ${{in_dir}}/*.bed`)

head -n 1 ${{files[1]}} > {output.table}

for perID_file in ${{files[@]}}
do
	sed -e '1d' $perID_file >> temp_out_table #{output.table}
done

#get intersection
bedtools intersect -a {input.rgn_bed} -b temp_out_table > temp_region_alignment

###TEMPORARY: NEED TO CHANGE TO ACCOUNT FOR BED FILES LACKING NAME COLUMN
cut -f 4 temp_region_alignment | paste - temp_out_table >> {output.table}
#### TEMPORARY

rm temp_region_alignment
rm temp_out_table
"""
