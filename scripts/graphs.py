from collections import Counter
import pandas as pd
import json
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import argparse

def GenesAndCodons(aa_muts, nt_muts):
    codons=[]
    all_genes =[] 
    with open(aa_muts) as a_muts:
        with open(nt_muts) as n_muts:
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
    genes =['F', 'G', 'M', 'M2', 'NS1', 'NS2', 'P', 'SH', 'N']
    muts_in_genes = dict()
    muts_in_genes_correct_index=dict()
    df= GenesAndCodons(aamutations, ntmutations)
    for gene in genes:
        muts_in_genes[gene]= df.loc[df['Gene']==gene]

    for gene, muts in muts_in_genes.items():
        muts = muts.reset_index(drop=True)
        muts_in_genes_correct_index[gene]=muts
    return(muts_in_genes_correct_index)


def AA_Mutations(aamutations, ntmutations):
    aa_m = dict()
    with open(aamutations) as f:
        with open(ntmutations) as g:
            genes = ('F','G','M','M2','NS1','NS2','P','SH','N')
            aamuts = json.load(f)
            for gene in genes:
                mut_list=[]
                for k, n in aamuts['nodes'].items():  
                            for i,j in n['aa_muts'].items():
                                if j!=[] and i ==gene:
                                    mut_list.append(j)
                flatlist =[item for sublist in mut_list for item in sublist]
                flatlist = [int(i[1:-1]) for i in flatlist]
                aa_m[gene]=flatlist
    return(aa_m)


def non_synonymous_or_synonymous(aa_muts, nt_muts):
    aa_mutations = AA_Mutations(aa_muts, nt_muts)
    mutations_in_genes = MutationsineachGene(aa_muts, nt_muts)
    synonymousmutations =[]
    nonsynonymousmutations =[]
    ratios=[]
    sel =[]
    listofgenes =('F','G','M','M2','NS1','NS2','P','SH','N')
    for gene in listofgenes:
        all_nonsyn_muts =[]
        all_syn_muts =[]
        for (gene_,mutation), (gene__,aa_mut) in zip(mutations_in_genes.items(), aa_mutations.items()):
            if gene_ == gene and gene__ == gene:
                a =list(mutation['Codon'])
                all_muts = Counter(a)
                amino_acid_muts = Counter(aa_mut)
                synonymous_muts = all_muts-amino_acid_muts
                

        for genes, j in amino_acid_muts.items():
            all_nonsyn_muts.append(j)
        nonsynonymousmutations.append(sum(all_nonsyn_muts))
        for k, l in synonymous_muts.items():
            all_syn_muts.append(l)
        synonymousmutations.append(sum(all_syn_muts))
        for a, b in zip(nonsynonymousmutations, synonymousmutations):
            ratio = a/b #ratio of nonsynonymous to synonymous mutations
            if ratio>1:selection =('adaptive')
            elif ratio<1: selection = ('purifying')
            elif ratio ==1: selection =('neutral')
        sel.append(selection)
        ratios.append(ratio)

    df = pd.DataFrame({"gene":listofgenes, "synonymous mutations": synonymousmutations, "nonsynonymous mutations":nonsynonymousmutations, "dN/dS ratio":ratios, "selection":sel })
    return(df)

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="analyse synonymous and nonsynonymous",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--aa', required=True, type=str, help="aa json file")
    parser.add_argument('--nt', type=str, help="nt json file")
    parser.add_argument('--output', type=str,  help="output graph png")
    parser.add_argument('--table', type=str,  help="output table csv")
    args = parser.parse_args()

    """ratio of nonsynonymous mutations in G to nonsynonymous in F is higher than synonymous G to synonymous F"""
    df1 = non_synonymous_or_synonymous(args.aa, args.nt)

    with open(args.aa) as f:
        gene_length=[]
        dictofgenes=dict()
        aamuts = json.load(f)
        keys = ('F','G','M','M2','NS1','NS2','P','SH','N')

        for gene in keys:
            for key, node in aamuts['annotations'].items():
                if key == gene:
                    loc_list=[]
                    loc_list = list(range(node['start'],node['end']+1))
            dictofgenes[gene]=loc_list
        for gene, loc in dictofgenes.items():
            gene_length.append(len(loc))
        df1['length of gene'] =gene_length
        df1['synonymous mutation/gene'] = df1['synonymous mutations']/df1['length of gene']
        df1['nonsynonymous mutation/gene']=df1['nonsynonymous mutations']/df1['length of gene']

    #ax1 = df1.plot.scatter(x='length of gene',y='synonymous mutations', title="RSV-A Synonymous mutations")
    #fig = ax1.get_figure()

    #for label, c in zip(df1['gene'], colorlist):
        #df1.plot.scatter(x='length of gene',y='synonymous mutations', ax=ax, s=50, linewidth=0.1, label=label, color=c)
    #fig.savefig("synonymous mutations RSV-A.png")


    plt.figure(figsize=(8,6))
    sp_names = df1['gene'].to_list()
    print(sp_names)
    colors_ = np.array(["green","blue","yellow","pink","black","orange","gray","cyan","magenta"])
    scatter = plt.scatter(df1['length of gene'], 
                df1['synonymous mutations'],
                s=150, c=colors_)

    plt.xlabel("Gene Length", size=20)
    plt.ylabel("synonymous mutations", size=20)
    plt.legend(handles=scatter.legend_elements()[0], 
            labels=sp_names,
            title="gene")
    plt.savefig(args.output)


    csv_file = non_synonymous_or_synonymous(args.aa, args.nt)
    csv_file.to_csv(args.table)