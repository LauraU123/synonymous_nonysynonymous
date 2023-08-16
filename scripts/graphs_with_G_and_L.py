from collections import Counter, defaultdict
import pandas as pd
import json, math, argparse
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO, Phylo


def GenesAndCodons(aa_muts, nt_muts):
    codons, all_genes= ([] for i in range(2))
    with open(aa_muts) as a_muts:
        with open(nt_muts) as n_muts:
            dictofgenes=dict()
            aamuts = json.load(a_muts)
            ntmuts = json.load(n_muts)
            genes = ('F','G','M','M2','NS1','NS2','P','SH','N', "L") #genes in RSV

            for gene in genes:
                for key, node in aamuts['annotations'].items():
                    if key == gene: location_of_gene = list(range(node['start'],node['end']+1)) #where each gene starts and ends
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
    genes =['F', 'G', 'M', 'M2', 'NS1', 'NS2', 'P', 'SH', 'N', "L"]
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
            genes = ('F','G','M','M2','NS1','NS2','P','SH','N', "L")
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

def scaled_mutations(reffile):
    read_ref = SeqIO.read(reffile, "genbank")
    cds_, sequence_ref_cds = (dict() for i in range(2))
    for feature in read_ref.features:
        if feature.type == "CDS": cds_[feature.qualifiers["gene"][0]] = list(feature.location)

    for gene, cds in cds_.items(): sequence_ref_cds[gene] = read_ref.seq[cds[0]:cds[-1]]
    translations = {'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'], 'E': ['GAA', 'GAG'], 'D': ['GAT', 'GAC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'N': ['AAT', 'AAC'], 'M': ['ATG'], 'K': ['AAA', 'AAG'], 'Y': ['TAT', 'TAC'], 'I': ['ATT', 'ATC', 'ATA'], 'Q': ['CAA', 'CAG'], 'F': ['TTT', 'TTC'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], '*': ['TAA', 'TAG', 'TGA'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'], 'H': ['CAT', 'CAC']}
    
    synonymous_possibilities, nonsynonymous_possibilities = (0 for i in range(2))
    for gene, sequence in sequence_ref_cds.items():
        for i, letter in enumerate(sequence):
            if i%3 == 0:
                codon = sequence[i: i+3]
                for  entry in translations.values():
                    if codon in entry:
                        synonymous = len(entry)
                        synonymous_possibilities += synonymous-1 
                        nonsynonymous = 9 - len(entry)-1
                        nonsynonymous_possibilities += nonsynonymous

    nonsyn_ratio = nonsynonymous_possibilities/(nonsynonymous_possibilities+synonymous_possibilities)
    syn_ratio = synonymous_possibilities/(nonsynonymous_possibilities+synonymous_possibilities)
    return(nonsyn_ratio, syn_ratio)

def non_synonymous_or_synonymous(reffile, aa_muts, nt_muts):
    non_ratio, syn_ratio = scaled_mutations(reffile)
    aa_mutations = AA_Mutations(aa_muts, nt_muts)
    mutations_in_genes = MutationsineachGene(aa_muts, nt_muts)
    synonymousmutations, nonsynonymousmutations, ratios, sel =([] for i in range(4))

    listofgenes =('F','G','M','M2','NS1','NS2','P','SH','N', "L")
    for gene in listofgenes:
        all_nonsyn_muts, all_syn_muts = ([] for i in range(2))
        for (gene_,mutation), (gene__,aa_mut) in zip(mutations_in_genes.items(), aa_mutations.items()):
            if gene_ == gene and gene__ == gene:
                a =list(mutation['Codon'])
                all_muts = Counter(a)
                amino_acid_muts = Counter(aa_mut)
                synonymous_muts = all_muts-amino_acid_muts
                
        for genes, j in amino_acid_muts.items(): all_nonsyn_muts.append(j)
        nonsynonymousmutations.append(sum(all_nonsyn_muts))
        for k, l in synonymous_muts.items(): all_syn_muts.append(l)
        synonymousmutations.append(sum(all_syn_muts))
        for a, b in zip(nonsynonymousmutations, synonymousmutations):
            ratio = (a/non_ratio)/(b/syn_ratio) #ratio of nonsynonymous to synonymous mutations
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
    parser.add_argument('--outputnonsyn', type=str, help="output file for nonsynonymous mutations")
    parser.add_argument('--ref', required=True, help="genbank reference file" )
    parser.add_argument('--table', type=str,  help="output table csv")
    parser.add_argument('--tree', type=str,  help="output table csv")
    parser.add_argument('--G_table', required=True, help="output G gene table file")
    args = parser.parse_args()

    tree_file  = Phylo.read(args.tree, "newick")
    tree_file.root_at_midpoint()
    tree_file.find_clades()

    #length of tree with and without duplication
    total_len = tree_file.total_branch_length()
    print(total_len)



    """ratio of nonsynonymous mutations in G to nonsynonymous in F is higher than synonymous G to synonymous F"""
    df1 = non_synonymous_or_synonymous(args.ref, args.aa, args.nt)
    df1.to_csv(args.table)
    ratio_nonsyn = scaled_mutations(args.ref)[0]
    ratio_syn = scaled_mutations(args.ref)[1]

    with open(args.aa) as f:
        gene_length=[]
        dictofgenes=dict()
        aamuts = json.load(f)
        keys = ('F','G','M','M2','NS1','NS2','P','SH','N', "L")

        for gene in keys:
            for key, node in aamuts['annotations'].items():
                if key == gene:
                    loc_list=[]
                    loc_list = list(range(node['start'],node['end']+1))
            dictofgenes[gene]=loc_list
        for gene, loc in dictofgenes.items():
            gene_length.append(len(loc))
        df1['length of gene'] =gene_length
        df1['synonymous mutation/gene'] = (df1['synonymous mutations']*3/df1['length of gene'])/total_len #per codon
        df1['nonsynonymous mutation/gene']=(df1['nonsynonymous mutations']*3/df1['length of gene'])/total_len #per codon

    plt.figure(figsize=(8,6))
    gene_names = df1['gene'].to_list()
    gene_name = df1['gene']
    colors_ = np.array(["palegreen","violet","yellow","pink","lightskyblue","orange","gray","cyan","bisque", "skyblue"])
    scatter = plt.plot(df1['gene'], 
                df1['synonymous mutation/gene'], 'o', markersize=10)
    for i in range(0, len(df1['length of gene'])):
        plt.text(df1['length of gene'][i] - 12, df1['synonymous mutation/gene'][i], f'{gene_name[i]}')
    plt.xlabel("Gene", size=14)
    plt.xticks(size=14)
    plt.yticks(size=14)
    plt.ylabel("synonymous mutations per codon", size=14)
    plt.savefig(args.output)

    plt.figure(figsize=(8,6))
    scatter_1 = plt.plot(df1['gene'], 
                df1['nonsynonymous mutation/gene'],
                'o', markersize=10)
    for i in range(0, len(df1['length of gene'])):
        plt.text(df1['length of gene'][i] - 12, df1['nonsynonymous mutation/gene'][i], f'{gene_name[i]}')
    plt.xlabel("Gene", size=14)
    plt.ylabel("nonsynonymous mutations per codon", size=14)
    plt.xticks(size=14)
    plt.yticks(size=14)
    plt.savefig(args.outputnonsyn)

    mylist=[]
    a = AA_Mutations(args.aa, args.nt)
    nt_muts = MutationsineachGene(args.aa, args.nt)['G']
    ntmuts_counter =(Counter(list(nt_muts['Codon'])))
    for j in a['G']: mylist.append(j)

    aa_muts_counter = (Counter(mylist))
    synonymous = ntmuts_counter-aa_muts_counter
    synonymouslist =[]
    for i, j in synonymous.items():
        for n in range(j): synonymouslist.append(i)
        
    with open(args.aa) as f: 
        aamuts = json.load(f)
        for key, node in aamuts['annotations'].items(): 
            if key == 'G':  g = list(range(node['start'],node['end']+1))

        G_gene_length = (math.floor(len(g)/3)+1)
        N_term = [*range(1,36)]
        T_Membrane = [*range(36,67)]
        Mucinlike_I = [*range(67,164)]
        Centralconserved = [*range(164,185)] 
        Heparin_binding = [*range(185,198)]
        Mucinlike_II = [*range(198, G_gene_length+1)]
        
        dictofdomains = {'N-terminal':N_term,'Transmembrane':T_Membrane,'Mucin-like I':Mucinlike_I,'Central conserved domain':Centralconserved,'Heparin-binding domain':Heparin_binding,'Mucin-like II': Mucinlike_II}
        dictionary_synonymous, dictionary_nonsynonymous = (defaultdict(list) for i in range(2))
        for name_, location_ in dictofdomains.items():
            for mut in synonymouslist:
                if mut in location_: dictionary_synonymous[name_].append(mut)
            for aa_mut in a['G']:
                if aa_mut in location_:  dictionary_nonsynonymous[name_].append(aa_mut)
        regions, nonsyn, syn, ratios =([] for i in range(4))
        for (name__, nonsynonymous), (name___, synonymous) in zip(dictionary_nonsynonymous.items(), dictionary_synonymous.items()):
            if name__ == name___: regions.append(name__), nonsyn.append(len(nonsynonymous)), syn.append(len(synonymous))
        for nonsynonymous_mutations, synonymous_mutations in zip(nonsyn, syn):
            ratio = (nonsynonymous_mutations/ratio_nonsyn)/(synonymous_mutations/ratio_syn)
            ratios.append(ratio)
        df = pd.DataFrame({'Region': regions,'Nonsynonymous Mutations': nonsyn, 'Synonymous Mutations':syn, "dN/dS ratio": ratios })
        df.to_csv(args.G_table)