params.genomes = ""
params.output_folder = ""
params.diamond_db = ""
params.seeds = ""
params.gvdb_tsv = ""
params.taxa_colors = ""

include { annotate; diamond_gvogs; get_markers; prepare_backbone; 
        concat_seeds_markers; run_mafft; run_divvier; run_trimal; 
        run_iqtree; color_tree; run_mafft_add} from 'gv-phylo/processes.nf'


workflow {
    genomes = Channel.fromPath(params.genomes)
    diamond_db = Channel.fromPath(params.diamond_db).first()
    seeds = Channel.fromPath(params.seeds)
    gvdb_tsv = Channel.fromPath(params.gvdb_tsv).first()
    taxa_colors = Channel.fromPath(params.taxa_colors).first()
    
    annotate(genomes)
    diamond_gvogs(annotate.out.proteomes, diamond_db)
    get_markers(diamond_gvogs.out.gvogs_out, seeds, annotate.out.proteomes)
    prepare_backbone(seeds, gvdb_tsv)
    marker_add = get_markers.out.markers.collectFile() { it ->
        [ "${it.name.tokenize('.')[1]}.additional.faa", it.text + '\n' ]
    }

    if (!params.aligned_seeds) {
      concat_seeds_markers(marker_add, prepare_backbone.out.backbone)
      run_mafft(concat_seeds_markers.out.concat_markers)
      run_divvier(run_mafft.out.alignments)
      run_trimal(run_divvier.out.divvied)
      run_iqtree(run_trimal.out.trimaled)
    else {
      aligned_seeds = Channel.fromPath(params.aligned_seeds)
      run_mafft_add(marker_add, aligned_seeds)
      run_iqtree(run_mafft_add.out.alignments)
    }
    color_tree(run_iqtree.out.treefiles, taxa_colors)
}
