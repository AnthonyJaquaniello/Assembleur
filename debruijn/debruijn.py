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
            d[kmer] = seq.count(kmer)
    return d

kmer_dict = build_kmer_dict(args.i, int(args.k))
print(kmer_dict)

## Construction de l'arbre de de Bruijn

def build_graph(dico):
    graph = nt.DiGraph()
    for key in dico:
        prefix = key[0:-1]
        suffix = key[1:]   
        graph.add_edge(prefix, suffix, weight = dico[key])
    return graph
l
graph = build_graph(kmer_dict)

    

