#! /usr/bin/env python3

from Bio import SeqIO
import re
import click

@click.command()
@click.argument('fasta', type=str, nargs=-1)
@click.option('--gff', type=str)

def main(fasta, gff):
    old2new = fix_gff(gff)
    for f in fasta:
        fix_fasta(f, old2new)

def fix_fasta(infile, old2new):
    recs = []
    for rec in SeqIO.parse(infile, 'fasta'):
        for o, n in old2new.items():
            if o in rec.description:
                rec.id = n.replace('ID=', '').replace(';', '')
        recs.append(rec)
    with open(infile, 'w') as out:
            SeqIO.write(recs, out, 'fasta')

def fix_gff(infile):
    lines = []
    old2new = {}
    for line in open(infile):
        line = line.strip()
        if not line.startswith("#"):
            llist = line.split()
            old_id = re.search('ID=([0-9]+)_([0-9]+);', line)
            new_id = re.sub(old_id.group(), "ID={}_{:03d};".format(llist[0], int(old_id.group(2))), old_id.group())
            line = re.sub(old_id.group(), new_id, line)
            line += "locus_tag={};".format(new_id)
            lines.append(line)
            old2new[old_id.group()] = new_id
    with open(infile, 'w') as out:
        for line in lines:
            print(line, file=out)
    return old2new
            
if __name__ == '__main__':
    main()
    
