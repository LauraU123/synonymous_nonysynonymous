def GenesAndCodons(aa_muts, nt_muts):

    import json
    import math
    import pandas as pd

    with open(aa_muts) as a_muts:
        with open(nt_muts) as n_muts:
            codons=[]
            indices=[] 
            dictofgenes=dict()
            aamuts = json.load(a_muts)
            ntmuts = json.load(n_muts)
            genes = ('F','G','M','M2','NS1','NS2','P','SH','N') #genes in RSV

            for gene in genes:
                for key, node in aamuts['annotations'].items():
                    if key == gene:
                        g=[]
                        g = list(range(node['start'],node['end']+1)) #where each gene starts and ends
                dictofgenes[gene]=g
        
            for k, n in aamuts['nodes'].items():
                for key, node in ntmuts['nodes'].items():
                    if k == key:
                        for y in node['muts']:
                            numbers =[]
                            number =int(y[1:-1])
                            numbers.append(number)
                            for x in numbers:
                                for a,g in dictofgenes.items():
                                    if x in g:
                                        index = a
                                        codon = (math.floor((x-g[0])/3))+1        
                                        indices.append(index)
                                        codons.append(codon)
        df=pd.DataFrame({'Gene':indices,'Codon':codons})
        return(df)