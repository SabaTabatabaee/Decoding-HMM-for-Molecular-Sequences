#! /usr/bin/env python

from HMM import *
import argparse
import sys

def read_FASTA(stream):
    '''
    This function reads a FASTA file from a given stream and returns a dictionary mapping identifiers to sequences
    '''
    seqs = {}; name = None; seq = ''
    for line in stream:
        l = line.strip()
        if len(l) == 0:
            continue
        if l[0] == '>':
            if name is not None:
                assert len(seq) != 0, "Malformed FASTA"
                seqs[name] = seq
            name = l[1:]
            assert name not in seqs, "Duplicate sequence ID: %s" % name
            seq = ''
        else:
            seq += l
    assert name is not None and len(seq) != 0, "Malformed FASTA"
    seqs[name] = seq
    return seqs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-m', '--model', required=True, type=str, help="HMM model")
    parser.add_argument('-q', '--query', required=True, type=str, help="Query sequences")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Write Viterbi score to this file")

    args = parser.parse_args()

    myHMM = HMM()
    myHMM.load(args.model)
    queries = read_FASTA(open(args.query,'r'))

    with open(args.output,'w') as fout:
        for q in queries:
            seq = queries[q]
            Vscore, aln = myHMM.Viterbi(seq)
            Vscore = round(Vscore, 5)
            fout.write(q + " " + str(Vscore) + " " + aln + "\n")
            #sys.exit()
