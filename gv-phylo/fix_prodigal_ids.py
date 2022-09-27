#! /home/mschoen/anaconda3/bin/python3

from Bio import SeqIO
import re
import click

@click.command()
@click.option('--faa', type=str)
@click.option('--faa_out', type=str)
@click.option('--ffn', type=str)
@click.option('--ffn_out', type=str)
@click.option('--gff', type=str)
@click.option('--gff_out', type=str)
@click.option('-p', '--prefix', type=str)

def main(faa, faa_out, ffn, ffn_out, gff, gff_out, prefix):
    fix_fasta(faa, faa_out, prefix)
    fix_fasta(ffn, ffn_out, prefix)
    fix_gff(gff, gff_out, prefix)
    
def fix_fasta(infile, outfile, prefix):
    with open(outfile, 'w') as out:
        for rec in SeqIO.parse(infile, 'fasta'):
            rec.description = re.sub("{} ".format(rec.id), "", rec.description)
            new_id = "{}_{}".format(prefix, rec.id)
            rec.description = re.sub('ID=[0-9_]+;', "ID={};".format(new_id), rec.description)
            rec.description += ";locus_tag={}".format(new_id)
            rec.id = new_id
            SeqIO.write(rec, out, 'fasta')

def fix_gff(infile, outfile, prefix):
    with open(outfile, 'w') as out:
        for line in open(infile):
            line = line.strip()
            if not line.startswith("#"):
                llist = line.split()
                line = re.sub(llist[0], "{}_{}".format(prefix, llist[0]), line)
                line = re.sub('ID=[0-9]+(_[0-9]+);', r"ID={}_{}\1;".format(prefix, llist[0]), line)
                new_id = re.search('ID=([A-Za-z0-9_]+);', line)
                line += "locus_tag={};".format(new_id.group(1))
            print(line, file=out)
            
if __name__ == '__main__':
    main()
    
