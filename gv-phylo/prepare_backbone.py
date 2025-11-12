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
@click.option('-s', '--subgroup', type=str, required=False)
@click.option('-i', '--ingroup', type=str, required=False)
@click.option('-e', '--exclude', type=str, required=False)
@click.option('--outgroup_per_family', type=int, default=1, required=False)
@click.option('--ingroup_per_family', type=int, default=5, required=False)


def main(gvdb, seeds, selectors, subgroup, ingroup, exclude, outgroup_per_family, ingroup_per_family):
    df = pd.read_csv(gvdb)
    df['Family'] = df.apply(lambda x: x['Family'] if x['Family_AKA'] == '-' else x['Family_AKA'], axis=1)
    df['Taxonomy'] = df.apply(lambda x: "_".join(x[['Class', 'Order', 'Family', 'clean_name']]), axis=1)
    df['Sequencing-approach'] = pd.Categorical(df['Sequencing-approach'], ["isolate", "cultivation-independent"], ordered=True)
    if subgroup:
        df = df[df['Taxonomy'].str.contains(subgroup)]
    selection = make_selection(df, selectors, ingroup, outgroup_per_family, ingroup_per_family)
    if exclude:
        for e in exclude.split(','):
            selection = selection[-selection['Taxonomy'].str.contains(e)]
    write_selection(selection, seeds)
    
def make_selection(df, selectors, ingroup, outgroup_per_family, ingroup_per_family):
    # select all with under 4 contigs:
    # sel = df[df['num_seqs'] < 4]
    # select 'best' 5 genomes per family, sorted by count and isolates first
    sel = pd.DataFrame(columns=df.columns)
    for d in df[df['Taxonomy'].str.contains(ingroup)].groupby('Family'):
        sel = pd.concat([sel,d[1].sort_values(by=['num_seqs', 'Sequencing-approach']).iloc[0:ingroup_per_family]])
    # get only one rep per family in outgroup
    for d in df[~df['Taxonomy'].str.contains(ingroup)].groupby('Family'):
        sel = pd.concat([sel,d[1].sort_values(by=['num_seqs', 'Sequencing-approach']).iloc[0:outgroup_per_family]])
    # manually add based on selection labels 
    for s in selectors:
        if ':' in s:
            s, num = s.split(':')
            sel = pd.concat([sel,df[df['Taxonomy'].str.contains(s)].sort_values(by=['num_seqs', 'Sequencing-approach']).iloc[0:int(num)]])
        else:
            sel = pd.concat([sel,df[df['Taxonomy'].str.contains(s)]])
    sel.drop_duplicates(inplace=True)
    return sel

def write_selection(sel, seeds):
    with open('{}.selection.faa'.format(os.path.basename(seeds).split('.')[0]), 'w') as out:
        for rec in SeqIO.parse(seeds, 'fasta'):
            genome = '_'.join(rec.id.split('..')[0].split('_')[0:-1])
            if genome in list(sel['clean_name']):
                prot = f"{rec.id.split('..')[0].split('_')[-1]}_{rec.id.split('..')[1]}"
                genome_tax = sel[sel['clean_name'] == genome]['Taxonomy'].item()
                rec.id = f'{genome_tax}..{prot}'
                rec.description = ""
                SeqIO.write(rec, out, 'fasta')

if __name__ == '__main__':
    main()
