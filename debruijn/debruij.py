import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i")
parser.add_argument("-k")
parser.add_argument("-o")
args = parser.parse_args()

def read_fastq(path):
    with open(path, "r") as filin:
        for line in filin:
            if line[0] in ["A","T","G","C"]:
                yield line[:-1]

def cut_kmer(seq, klen):
    for i in range(0, len(seq)-klen+1):
        yield seq[i:i+klen]

def build_kmer_dict(i, k):
    list_d = []
    seq_gen = read_fastq(i)
    for seq in seq_gen:
        d = {}
        for kmer in cut_kmer(seq, k):
            d[kmer] = seq.count(kmer)
        list_d.append(d)
    return list_d

liste_dico = build_kmer_dict(args.i, int(args.k))



    

