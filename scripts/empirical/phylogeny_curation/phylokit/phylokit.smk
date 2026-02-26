# Authors: Kyle Gontjes & Dhatri Badri
# Date: 02-20-2024

configfile: "config/config.yaml"

import pandas as pd
import os
import re
 
PREFIX = config["prefix"]

print(f"Config: {config}") 
print(f"PREFIX: {PREFIX}")

if not os.path.exists("results/"):
    try:
        os.makedirs("results/")
    except OSError as e:
        print(f"Error creating directory: {e}")

# Rule all information
rule all:
    input:
        recombination_predictions_gff = expand("results/gubbins/{prefix}.recombination_predictions.gff",prefix=PREFIX), 
        gubbins_masked = expand("results/gubbins_masked/{prefix}_gubbins_masked.fa",prefix=PREFIX),
        masked_recomb_pos = expand("results/gubbins_masked/{prefix}_masked_recomb_positions.txt",prefix=PREFIX),
        gubbins_masked_var_sites = expand("results/gubbins_masked/{prefix}_gubbins_masked_var_sites.fa",prefix=PREFIX), 
        ml_tree_newick_fmt = expand("results/IQtree/{prefix}.treefile",prefix=PREFIX)

rule gubbins:
    input: 
        alignment = config["alignment"]
    output:
        recombination_predictions_gff = f"results/gubbins/{{prefix}}.recombination_predictions.gff"
    params:
        pwd = config["gubbins_dir"],
        file_prefix = config["prefix"],
        threads = config["threads"],
        first_tree_algorithm = config["first_tree_algorithm"],
        first_model = config["first_model"], 
        outgroup = config["outgroup"] 
    singularity:
        "docker://staphb/gubbins:3.3.5"
    shell:
        """
        cd {params.pwd} && \
        run_gubbins.py --prefix {params.file_prefix} --threads {params.threads} --verbose --no-cleanup --first-tree-builder {params.first_tree_algorithm} --first-model {params.first_model} --best-model --outgroup {params.outgroup} {input.alignment}
        """

rule mask_gubbins:
    input:
        recombination_predictions_gff = lambda wildcards: expand(f"results/gubbins/{wildcards.prefix}.recombination_predictions.gff")
    output:
        gubbins_masked = f"results/gubbins_masked/{{prefix}}_gubbins_masked.fa",
        masked_recomb_pos = f"results/gubbins_masked/{{prefix}}_masked_recomb_positions.txt",
        gubbins_masked_var_sites = f"results/gubbins_masked/{{prefix}}_gubbins_masked_var_sites.fa"
    params:
        outdir = "results/gubbins_masked",
        prefix =  config["prefix"],
        alignment = config["alignment"]
    wrapper:
        "file:wrapper_functions/mask_gubbins" 

rule iqtree:
    input:
        gubbins_masked_var_sites = lambda wildcards: expand(f"results/gubbins_masked/{wildcards.prefix}_gubbins_masked_var_sites.fa") 
    output:
        ml_tree_newick_fmt = f"results/IQtree/{{prefix}}.treefile" 
    params:
       iqtree_model = config["iqtree_model"],
       bootstrap_count = config["bootstrap_count"],
       outgroup = config["outgroup"],
       num_unsuccessful_iterations = config["num_unsuccessful_iterations"],
       iqtree_prefix = expand(str(config["iqtree_dir"]) + "/" + str(config["prefix"]))
    conda:
        "envs/iqtree.yaml"
    log:
        iqtree_log = "logs/IQtree/{prefix}.iqtree.log"
    shell:
        "iqtree -s {input.gubbins_masked_var_sites} -nt AUTO -m {params.iqtree_model} -B {params.bootstrap_count} -o {params.outgroup} -nstop {params.num_unsuccessful_iterations} -pre {params.iqtree_prefix} &> {log.iqtree_log}" 
