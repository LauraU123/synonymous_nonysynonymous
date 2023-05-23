
SUBTYPES = ["a", "b"]
rule all:
    input:
        expand("results/{subtype}/synonymous.png", subtype=SUBTYPES)

rule graphs:
    input:
        aa = "data/{subtype}/aa_muts.json",
        nt = "data/{subtype}/nt_muts.json"
    output:
        graph = "results/{subtype}/synonymous.png",
        graphnonsyn = "results/{subtype}/nonsynonymous.png",
        table = "results/{subtype}/muts.csv"
    shell:
        """
        python3 scripts/graphs.py \
        --aa {input.aa} \
        --nt {input.nt} \
        --output {output.graph} \
        --outputnonsyn {output.graphnonsyn} \
        --table {output.table}
        """
