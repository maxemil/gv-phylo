#! /usr/bin/env python3

from Bio import SeqIO
import click
import os
import numpy as np

@click.command()
@click.option('--alignment', '-a', type=str, required=True)
@click.option('--output', '-o', type=str, required=True)
@click.option('--fraction', '-f', type=float, required=True)


def main(alignment, output, fraction):
    seqdict = SeqIO.to_dict(SeqIO.parse(alignment, 'fasta'))
    median_len = np.median([len(ungap(s)) for s in seqdict.values()])
    with open(output, 'w') as out:
        for rec in seqdict.values():
            if (len(ungap(rec)) / median_len) >= fraction:
                SeqIO.write(rec, out, 'fasta')
            else:
                print(f"removed {rec.id}, as it was less than {fraction} of the median sequence length.")

def ungap(rec):
    # undef_chars = set(['-', 'X'])
    return rec.seq.replace('-', '').replace('X', '')

if __name__ == '__main__':
        main()

