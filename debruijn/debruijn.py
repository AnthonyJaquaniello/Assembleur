import os
import argparse
import networkx as nt
from statistics import stdev, mean
from networkx.algorithms.shortest_paths.generic import shortest_path

def select_best_path():
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
            if line[0] in ["A", "T", "G", "C"]:
                yield line[:-1]

def cut_kmer(seq, klen):
    for i in range(0, len(seq)-klen+1):
        yield seq[i: i+klen]

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

kmer_dict = build_kmer_dict("data/eva71_two_reads.fq", int(21))
#../data/eva71_two_reads.fq
## Construction de l'arbre de de Bruijn

def build_graph(dico):
    graphi = nt.DiGraph()
    for key in dico:
        prefix = key[:-1]
        suffix = key[1:]
        graphi.add_edge(prefix, suffix, weight=dico[key])
    return graphi

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

output_node = get_sink_nodes(graph)

def get_contigs(graph, inode, outnode):
    tuple_list = []
    for nodin in inode:
        for nodout in outnode:
            if len(list(nt.all_simple_paths(graph, nodin, nodout))) != 0:
                path = shortest_path(graph, nodin, nodout)
                contig = path[0]
                for i in range(1, len(path)):
                    contig += path[i][-1]
                tuples = []
                tuples.append(contig)
                tuples.append(len(contig))
                tuple_list.append(tuple(tuples))
    return tuple_list

#contig_list = get_contigs(graph, input_node, output_node)

def fill(text, width=80):
    return os.linesep.join(text[i: i+width] for i in range(0, len(text), width))

def save_contigs(tuple_list, path):
    with open(path, "w") as filin:
        for i in range(len(tuple_list)):
            filin.write(">contig_{} len={}\n".format(i, tuple_list[i][1]))
            filin.write(fill(tuple_list[i][0]))
            filin.write("\n")

#save_contigs(contig_list, "save_contig.txt")

def std(liste):
    return stdev(liste)

path_gen = list(nt.all_simple_paths(graph, input_node[0], output_node[0]))

def path_average_weight(graph, path):
    path = list(path)
    all_w = []
    for i in range(0, len(path)-1):
        for key in graph:
            if key == path[i]:
                try:
                    d = dict(graph[key])
                    w = list(d.values())[0]["weight"]
                    all_w.append(w)
                except:
                    all_w.append(0)
    return mean(all_w)

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if delete_entry_node:
            graph.remove_node(path[0])
        if delete_sink_node:
            graph.remove_node(path[-1])
        for i in range(1, len(path)-1):
            graph.remove_node(path[i])
    return graph

def select_best_path(graph, path_list, len_list, weight_list, delete_entry_node, delete_sink_node):
    best_path = path_list[index(max(weight_list))]
    path_list.remove(path_list)
    remove_paths(path_list)
    pass 
