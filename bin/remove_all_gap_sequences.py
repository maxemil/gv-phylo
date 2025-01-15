#! /usr/bin/env python3

from Bio import SeqIO
import click
import os

@click.command()
@click.option('--alignment', '-a', type=str, required=True)
@click.option('--output', '-o', type=str, required=True)


def main(alignment, output):
    undef_chars = set(['-', 'X'])
    with open(output, 'w') as out:
        for rec in SeqIO.parse(alignment, 'fasta'):
            if not set(rec.seq).issubset(set(['-', 'X'])):
                SeqIO.write(rec, out, 'fasta')
            else:
                print("removed {}, as it was only made up of - and X.".format(rec.id))

if __name__ == '__main__':
        main()

