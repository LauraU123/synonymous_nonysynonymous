def GenesAndCodons(aa_muts, nt_muts):

    import json
    import math
    import pandas as pd

    with open(aa_muts) as a_muts:
        with open(nt_muts) as n_muts:
            codons=[]
            all_genes =[] 
            dictofgenes=dict()
            aamuts = json.load(a_muts)
            ntmuts = json.load(n_muts)
            genes = ('F','G','M','M2','NS1','NS2','P','SH','N') #genes in RSV

            for gene in genes:
                for key, node in aamuts['annotations'].items():
                    if key == gene:
                        location_of_gene=[]
                        location_of_gene = list(range(node['start'],node['end']+1)) #where each gene starts and ends
                dictofgenes[gene]=location_of_gene
        
            for k, n in aamuts['nodes'].items():
                for key, node in ntmuts['nodes'].items():
                    if k == key:
                        for y in node['muts']:
                            numbers =[]
                            number =int(y[1:-1])
                            numbers.append(number)
                            for pos in numbers:
                                for gene, location_of_gene in dictofgenes.items():
                                    if pos in location_of_gene:
                                        codon = (math.floor((pos-location_of_gene[0])/3))+1        
                                        all_genes.append(gene)
                                        codons.append(codon)
        df=pd.DataFrame({'Gene':all_genes,'Codon':codons})
        return(df)
    
def MutationsineachGene(aamutations, ntmutations):

    import pandas as pd
    genes =['F', 'G', 'M', 'M2', 'NS1', 'NS2', 'P', 'SH', 'N']
    dictionary = dict()
    dictionary1=dict()
    df=(GenesAndCodons(aamutations, ntmutations))
    for i in genes:

        dictionary[i]= df.loc[df['Gene']==i]

    for i,j in dictionary.items():

        a =(j['Codon'].value_counts())
        b = pd.DataFrame({'Codon':a.index, 'Frequency':a.values})

    for i,j in dictionary.items():
        j =j.reset_index(drop=True)
        dictionary1[i]=j
    return(dictionary1)


def AA_Mutations(aamutations, ntmutations):
    aa_m = dict()

    with open(aamutations) as f:
        with open(ntmutations) as g:

            lis1 =[]
            keys = ('F','G','M','M2','NS1','NS2','P','SH','N')
            import json
            aamuts = json.load(f)
            ntmuts = json.load(g)

            for a in keys:

                lis1=[]

                for k, n in aamuts['nodes'].items():  
                            for i,j in n['aa_muts'].items():
                                if j!=[] and i ==a: lis1.append(j)
                flatlist =[item for sublist in lis1 for item in sublist]
                #print(len(flatlist))
                flatlist = [int(i[1:-1]) for i in flatlist]
                aa_m[a]=flatlist
    return(aa_m)


def non_synonymous_or_synonymous(aa_muts, nt_muts):
    aa = AA_Mutations(aa_muts, nt_muts)
    mu = MutationsineachGene(aa_muts, nt_muts)
    from collections import Counter
    import pandas as pd
    synonymousmutations =[]
    nonsynonymousmutations =[]
    ratios=[]
    sel =[]
    listofgenes =('F','G','M','M2','NS1','NS2','P','SH','N')
    for gene in listofgenes:
        list1 =[]
        list2 =[]
        list3 =[]
        for (k,l), (i,j) in zip(mu.items(), aa.items()):
            if k == gene and i == gene:
                a =(list(l['Codon']))
                c = Counter(a)
                b = Counter(j)
                d = c-b
                

        for i, j in b.items():
            list1.append(j)
        nonsynonymousmutations.append(sum(list1))
        for k, l in d.items():
            list2.append(l)
        synonymousmutations.append(sum(list2))
        for a, b in zip(nonsynonymousmutations, synonymousmutations):
            ratio = a/b
            if ratio>1:selection =('adaptive')
            elif ratio<1: selection = ('purifying')
            elif ratio ==1: selection =('neutral')
        sel.append(selection)
        ratios.append(ratio)

    df = pd.DataFrame({"gene":listofgenes, "synonymous mutations": synonymousmutations, "nonsynonymous mutations":nonsynonymousmutations, "dN/dS ratio":ratios, "selection":sel })
    return(df)


"""ratio of nonsynonymous mutations in G to nonsynonymous in F is higher than synonymous G to synonymous F"""
df1 = non_synonymous_or_synonymous('results/b/genome/aa_muts.json', 'results/b/genome/nt_muts.json')
#print(df1['synonymous mutations'])

import json
with open('results/b/genome/aa_muts.json') as f:
    list1=[]
    dictofgenes=dict()
    aamuts = json.load(f)
    keys = ('F','G','M','M2','NS1','NS2','P','SH','N')

    for i in keys:
        for key, node in aamuts['annotations'].items():
            if key == i:
                g=[]
                g = list(range(node['start'],node['end']+1))
        dictofgenes[i]=g
    for (i, j) in dictofgenes.items():
        print(len(j))
        list1.append(len(j))
    df1['length of gene'] =list1
    df1['synonymous mutation/gene'] = df1['synonymous mutations']/df1['length of gene']
    df1['nonsynonymous mutation/gene']=df1['nonsynonymous mutations']/df1['length of gene']
    print(df1)

import matplotlib.pyplot as plt
from matplotlib import colors
#ax1 = df1.plot.scatter(x='length of gene',y='synonymous mutations', title="RSV-A Synonymous mutations")
#fig = ax1.get_figure()

#for label, c in zip(df1['gene'], colorlist):
    #df1.plot.scatter(x='length of gene',y='synonymous mutations', ax=ax, s=50, linewidth=0.1, label=label, color=c)
#fig.savefig("synonymous mutations RSV-A.png")


plt.figure(figsize=(8,6))
sp_names = df1['gene'].to_list()

scatter = plt.scatter(df1['length of gene'], 
            df1['synonymous mutations'],
            s=150)
plt.xlabel("Gene Length", size=20)
plt.ylabel("synonymous mutations", size=20)
# add legend to the plot with names
plt.legend(handles=scatter.legend_elements()[0], 
           labels=sp_names,
           title="gene")