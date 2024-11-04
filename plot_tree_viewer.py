#! /usr/bin/env python3

import click
import subprocess

@click.command()
@click.option('-t', '--treefile', type=str, required=True)
@click.option('-o', '--outprefix', type=str, required=True)
@click.option('-s', '--style', type=str, required=True, default='basic_tree')


def main(treefile, outprefix, style):
    if style == 'basic_tree':
        tree_viewer_commands = basic_tree(treefile, outprefix)
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

if __name__ == '__main__':
    main()