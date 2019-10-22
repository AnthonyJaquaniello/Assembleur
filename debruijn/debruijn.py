import argparse
import networkx as nt

def get_starting_nodes():
    pass

def std():
    pass

def get_sink_nodes():
    pass


def path_average_weight():
    pass


def remove_paths():
    pass


def select_best_path():
    pass


def save_contigs():
    pass


def get_contigs():
    pass


def solve_bubble():
    pass


def simplify_bubbles():
    pass


def solve_entry_tips():
    pass


def solve_out_tips():
    pass

parser = argparse.ArgumentParser()
parser.add_argument("-i")
parser.add_argument("-k")
parser.add_argument("-o")
args = parser.parse_args()

#Création du Graph de de Bruijn
## Identification des kmers uniques

def read_fastq(path):
    with open(path, "r") as filin:
        for line in filin:
            if line[0] in ["A","T","G","C"]:
                yield line[:-1]

def cut_kmer(seq, klen):
    for i in range(0, len(seq)-klen+1):
        yield seq[i:i+klen]

def build_kmer_dict(i, k):
    seq_gen = read_fastq(i)
    d = {}
    for seq in seq_gen:
        for kmer in cut_kmer(seq, k):
            if kmer not in d.keys():
                d[kmer] = 1
            else:
                d[kmer] += 1
    return d

kmer_dict = build_kmer_dict("data/eva71_two_reads.fq", 21)

## Construction de l'arbre de de Bruijn

def build_graph(dico):
    graph = nt.DiGraph()
    for key in dico:
        prefix = key[0:-1]
        suffix = key[1:]   
        graph.add_edge(prefix, suffix, weight = dico[key])
    return graph

graph = build_graph(kmer_dict)

## Parcours graph de De Bruijn
def get_starting_nodes(graph):
    l = []
    for d in graph.pred:
        if len(graph.pred[d]) == 0:
            l.append(d)
    return l

input_node = get_starting_nodes(graph)

def get_sink_nodes(graph):
    l = []
    for d in graph.succ:
        if len(graph.succ[d]) == 0:
            l.append(d)
    return l

print(input_node)
output_node = get_sink_nodes(graph)
print(output_node)    
