import os

workdir: "ortho_results" #add output directory

rule all:
    input:
        final = "results.txt"  #this is an input to rule all step

rule first_step:
    input:
        #nothing: can take inputs but for this we're not
        # rules don't have to have inputs
    output:
        txt = "a.txt", #names only need to be unique within rules: scope of name is within rule
        txt2 = "b.txt", #can always end input or output lists with comma and it won't screw up the interpreter. 
    shell:"""
echo "line 1" > {output.txt} #squigle braces enables me to generalize what I'm writing to for output. Uses python syntax to pass to shell command. use squiggles to access input or output
touch {output.txt2}
"""
rule result:
    input:
        txt = rules.first_step.output.txt #this passes output from rule first_step : equivalent to : <txt = "a.txt">
    output:
        result = "results.txt"
    run:
        #run python code
        f = open( input.txt )
        counter = 0
        for line in f:
            counter += 1
        open( output.result , "w+").write("{}\n".format(counter))
