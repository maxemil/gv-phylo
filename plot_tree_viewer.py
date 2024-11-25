#! /usr/bin/env python3

import click
import subprocess

@click.command()
@click.option('-t', '--treefile', type=str, required=True)
@click.option('-o', '--outprefix', type=str, required=True)
@click.option('-q', '--queryprefix', type=str, required=False)
@click.option('-s', '--style', type=str, required=False, default='basic_tree')


def main(treefile, outprefix, style, queryprefix):
    if style == 'basic_tree':
        tree_viewer_commands = basic_tree(treefile, outprefix)
    elif style == 'highlight_tree':
        tree_viewer_commands = highlight_tree(treefile, outprefix, queryprefix)
    else:
        print('couldnt read tree style argument, using basic_tree')
        tree_viewer_commands = basic_tree(treefile, outprefix)
    p = subprocess.run(['TreeViewerCommandLine'], input=tree_viewer_commands, capture_output=True, text=True)

def basic_tree(input, output_prefix):
    return f"""
open {input}
y

module enable Rooted tree style

module enable Reroot tree
y
option select #1
option set Mid-Point
update

module select #8
option select #1
option set Support
update

binary modules loaded {output_prefix}.tbi
y
y
pdf {output_prefix}.pdf
""".format(input, output_prefix)

def highlight_tree(input, output_prefix, highlight_prefix):
    return f"""
open {input}
y

module enable Rooted tree style

module enable Reroot tree
y
option select #1
option set Mid-Point
update

module select #8
option select #1
option set Support
update

module enable Replace attribute
option select #1
option set Name
option select #3
option set {highlight_prefix}
option select #8
option set Query
option select #9
option set Number
option select #10
option set 20
update

module enable Node shapes
option select Size
option set 0
option set attribute number Query
option select Auto fill colour by node
option set false
option select Stroke thickness
option set 0
update 


binary modules loaded {output_prefix}.tbi
y
y
pdf {output_prefix}.pdf
""".format(input, output_prefix, highlight_prefix)

if __name__ == '__main__':
    main()