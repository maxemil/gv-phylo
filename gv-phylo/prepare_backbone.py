#! /usr/bin/env python3

from Bio import SeqIO
import pandas as pd
import os
import re
import click 

@click.command()
@click.argument('gvdb', type=str, required=1)
@click.argument('seeds', type=str, required=1)
@click.argument('selectors', type=str, nargs=-1)
@click.option('-s', '--subgroup', type=str)
@click.option('-i', '--ingroup', type=str)


def main(gvdb, seeds, selectors, subgroup, ingroup):
    df = pd.read_csv(gvdb, sep='\t')
    df['Family'] = df.apply(lambda x: x['Family'] if x['Family_AKA'] == '-' else x['Family_AKA'], axis=1)
    df['Taxonomy'] = df.apply(lambda x: "_".join(x[['Class', 'Order', 'Family', 'common_name']]), axis=1)
    if subgroup:
        df = df[df['Taxonomy'].str.contains(subgroup)]
    selection = make_selection(df, selectors, ingroup)
    write_selection(selection, seeds)
    
def make_selection(df, selectors, ingroup):
    # select all with under 4 contigs:
    # sel = df[df['num_seqs'] < 4]
    # select 'best' 5 genomes per family, sorted by count and isolates first
    sel = pd.DataFrame(columns=df.columns)
    for d in df[df['Taxonomy'].str.contains(ingroup)].groupby('Family'):
        sel = pd.concat([sel,d[1].sort_values(by=['num_seqs', 'Sequencing-approach']).iloc[0:5]])
    # get only one rep per family in outgroup
    for d in df[~df['Taxonomy'].str.contains(ingroup)].groupby('Family'):
        sel = pd.concat([sel,d[1].sort_values(by=['num_seqs', 'Sequencing-approach']).iloc[0:1]])
    # manually add based on selection labels 
    for s in selectors:
         sel = pd.concat([sel,df[df['Taxonomy'].str.contains(s)]])
    sel.drop_duplicates(inplace=True)
    return sel

def write_selection(sel, seeds):
    p = re.compile('\.|\ |:|=|-')
    with open('{}.selection.faa'.format(os.path.basename(seeds).split('.')[0]), 'w') as out:
        for rec in SeqIO.parse(seeds, 'fasta'):
            genome = rec.id.split('|')[0].replace('.fna', '').replace('.fa', '').replace('.fltr', '')
            if genome in list(sel['genome_id']):
                if len(rec.id.split('|')) > 1:
                    prot = p.sub('_', rec.id.split('|')[1]).strip('_')
                else:
                    prot = 'unknown_protein'
                genome_tax = "_".join([genome, sel[sel['genome_id'] == genome]['Taxonomy'].item()])
                rec.id = p.sub('_', genome_tax)
                rec.id = f'{rec.id}..{prot}'
                rec.description = ""
                SeqIO.write(rec, out, 'fasta')

if __name__ == '__main__':
    main()
