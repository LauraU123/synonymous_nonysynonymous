
SUBTYPES = ["a", "b"]
rule all:
    input:
        expand("results/{subtype}/synonymous.png", subtype=SUBTYPES)

rule graphs:
    input:
        aa = "data/{subtype}/aa_muts.json",
        nt = "data/{subtype}/nt_muts.json"
    output:
        "results/{subtype}/synonymous.png"
    shell:
        """
        python3 scripts/graphs.py \
        --aa {input.aa} \
        --nt {input.nt} \
        --output {output}
        """