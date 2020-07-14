import os
SMK_DIR = os.path.dirname(workflow.snakefile) #workflow.snakefile == full path to snakefile 
shell.prefix("source {}/env.cfg; set -eo pipefail;".format(SMK_DIR)) #source env config before anything gets run ##these eo sets tells pipeline to stop running if one of lines fails.
configfile: "ortho.yaml"
SMS = list(config.keys()) #just strings of element names in yaml ["a","b"]
