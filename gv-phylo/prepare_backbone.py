#! /usr/bin/env python3

from Bio import SeqIO
import pandas as pd
import os
import click 

@click.command()
@click.argument('gvdb', type=str, required=1)
@click.argument('seeds', type=str, required=1)
@click.argument('selectors', type=str, nargs=-1)
@click.option('-s', '--subgroup', type=str)

def main(gvdb, seeds, selectors, subgroup):
    df = pd.read_csv(gvdb, sep='\t')
    df['Family'] = df.apply(lambda x: x['Family'] if x['Family_AKA'] == '-' else x['Family_AKA'], axis=1)
    df['Taxonomy'] = df.apply(lambda x: "_".join(x[['Class', 'Order', 'Family', 'common_name']]), axis=1)
    df = df[df['Taxonomy'].str.contains(subgroup)]
    selection = make_selection(df, selectors)
    write_selection(selection, seeds)
    
def make_selection(df, selectors):
    # select all with under 4 contigs:
    # sel = df[df['num_seqs'] < 4]
    # select 'best' 5 genomes per family, sorted by count and isolates first
    sel = pd.DataFrame(columns=df.columns)
    for d in df.groupby('Family'):
        sel = pd.concat([sel,d[1].sort_values(by=['num_seqs', 'Sequencing-approach']).iloc[0:5]])
    # manually add based on selection labels 
    for s in selectors:
         sel = pd.concat([sel,df[df['Taxonomy'].str.contains(s)]])
    sel.drop_duplicates(inplace=True)
    return sel

def write_selection(sel, seeds):
    with open('{}.selection.faa'.format(os.path.basename(seeds).split('.')[0]), 'w') as out:
        for rec in SeqIO.parse(seeds, 'fasta'):
            genome = rec.id.split('|')[0].replace('.fna', '').replace('.fa', '').replace('.fltr', '')
            if genome in list(sel['genome_id']):
                rec.id = "_".join([rec.id, sel[sel['genome_id'] == genome]['Taxonomy'].item()])
                rec.id = rec.id.replace(' ', '_').replace('|', '..').replace(':', '_').strip('_')
                rec.description = ""
                SeqIO.write(rec, out, 'fasta')

if __name__ == '__main__':
    main()
