
SUBTYPES = ["a", "b"]
rule all:
    input:
        expand("results/{subtype}/synonymous.png", subtype=SUBTYPES)

rule graphs:
    input:
        aa = "data/{subtype}/aa_muts.json",
        nt = "data/{subtype}/nt_muts.json",
        ref = "data/{subtype}/{subtype}reference.gbk",
        tree = "data/{subtype}/tree.nwk"
    output:
        graph = "results/{subtype}/synonymous.png",
        graphnonsyn = "results/{subtype}/nonsynonymous.png",
        table = "results/{subtype}/muts.csv",
        G = "results/{subtype}/G_muts.csv"
    shell:
        """
        python3 scripts/graphs_with_G_and_L.py \
        --aa {input.aa} \
        --nt {input.nt} \
        --tree {input.tree} \
        --output {output.graph} \
        --outputnonsyn {output.graphnonsyn} \
        --table {output.table} \
        --ref {input.ref} \
        --G_table {output.G}
        """
