params.genomes = ""
params.output_folder = ""
params.diamond_db = ""
params.seeds = ""
params.gvdb_tsv = ""
params.taxa_colors = ""
params.aligned_seeds = ""
params.treemode = 'iqtree-fast'

include { annotate; diamond_gvogs; get_markers; prepare_backbone; concat_seeds_markers; run_mafft; run_divvier; run_trimal; run_iqtree; color_tree; run_mafft_add} from './gv-phylo/processes.nf'

workflow {
    genomes = Channel.fromPath(params.genomes)
    diamond_db = Channel.fromPath(params.diamond_db).first()
    seeds = Channel.fromPath(params.seeds)
    gvdb_tsv = Channel.fromPath(params.gvdb_tsv).first()
    taxa_colors = Channel.fromPath(params.taxa_colors).first()
    
    annotate(genomes)
    diamond_gvogs(annotate.out.proteomes, diamond_db)
    get_markers(diamond_gvogs.out.gvogs_out, seeds, annotate.out.proteomes)
    marker_add = get_markers.out.markers.collectFile() { it ->
        [ "${it.name.tokenize('.')[1]}.additional.faa", it.text + '\n' ]
    }

    if (!params.aligned_seeds) {
      prepare_backbone(seeds, gvdb_tsv)
      concat_seeds_markers(marker_add, prepare_backbone.out.backbone)
      run_mafft(concat_seeds_markers.out.concat_markers)
      run_divvier(run_mafft.out.alignments)
      run_trimal(run_divvier.out.divvied)
      run_iqtree(run_trimal.out.trimaled)
    } else {
      aligned_seeds = Channel.fromPath(params.aligned_seeds)
      run_mafft_add(marker_add, aligned_seeds)
      run_iqtree(run_mafft_add.out.alignments)
    }
    color_tree(run_iqtree.out.treefiles, taxa_colors)
}

workflow.onComplete {
    File file = new File("$params.output_folder/${workflow.start}.log")
    file.append("Pipeline $workflow.scriptName started at $workflow.start \n")
    file.append("Pipeline $workflow.scriptName with hash $workflow.scriptId \n")
    file.append("Pipeline $workflow.scriptName was launched at $workflow.launchDir \n")
    file.append("Pipeline $workflow.scriptName was launched as $workflow.commandLine \n")
    file.append("Pipeline $workflow.scriptName completed at $workflow.complete \n")
    file.append("Execution status: ${ workflow.success ? 'OK' : 'failed' } \n")
}
