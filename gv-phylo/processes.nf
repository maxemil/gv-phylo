params.treemode = 'iqtree-fast'
params.marker_selection = "all"
params.selectors = 'Cafeteria'
params.ingroup = 'Megaviricetes'
params.subgroup = ''

process annotate {
  input:
    file genome
  output:
    path "${genome.simpleName}.faa", emit: proteomes
    path "${genome.simpleName}.{gff,ffn}", emit: prodigal_files
    
  publishDir "${params.output_folder}/annotation", mode: 'copy', pattern: '*.{gff,ffn,faa}'

  script:
    """
    prodigal -i $genome \
            -o ${genome.simpleName}.gff \
            -f gff \
            -d ${genome.simpleName}.ffn \
            -a ${genome.simpleName}.faa \
            -n -p meta
    fix_prodigal_ids.py \
            --gff ${genome.simpleName}.gff \
            ${genome.simpleName}.ffn \
            ${genome.simpleName}.faa
            
#    prokka --kingdom Viruses ${genome} \
#            --outdir annotation \
#            --prefix ${genome.simpleName} \
#            --locustag ${genome.simpleName}
    """
}

process diamond_gvogs {
  input:
    file proteome
    file diamond_db

  output:
    path "${proteome.simpleName}.${diamond_db.simpleName}.o6", emit: gvogs_out

  label "mid_cpu"
  publishDir "${params.output_folder}/diamond_gvogs", mode: 'copy'


  script:
    """
    diamond blastp --query ${proteome} \
          --db ${diamond_db} \
          --out ${proteome.simpleName}.${diamond_db.simpleName}.o6 \
          -f 6 -p ${task.cpus} 
    """
}

process get_markers {
  input:
    file outfile
    each seeds
    each proteome

  when:
    "${outfile.simpleName}" == "${proteome.simpleName}"

  output:
    path "${outfile.simpleName}.${seeds.simpleName}.fasta", emit: markers, optional: true

  publishDir "${params.output_folder}/markers", mode: 'copy'

  script:
    """
    #!/home/mschoen/anaconda3/bin/python3
    
    from Bio import SeqIO
    
    seeds = [rec.id for rec in SeqIO.parse('${seeds}', 'fasta')]
    
    prots = set()
    for line in open('${outfile}'):
        line = line.strip().split('\t')
        if line[1] in seeds:
            prots.add(line[0])
    if prots:
        with open("${outfile.simpleName}.${seeds.simpleName}.fasta", 'w') as out:
            if "${params.marker_selection}" == 'all':
                for rec in SeqIO.parse("${proteome}", 'fasta'):
                    if rec.id in prots:
                        SeqIO.write(rec, out, 'fasta')
            elif "${params.marker_selection}" == 'greedy':
                best_rec = None
                longest = 0
                for rec in SeqIO.parse("${proteome}", 'fasta'):
                    if rec.id in prots and len(rec.seq) > longest:
                        best_rec = rec
                        longest = len(rec.seq)
                SeqIO.write(best_rec, out, 'fasta')
    """
}

process prepare_backbone {
  input:
    file seeds
    file gvdg_genomes_tsv

  output:
    path "${seeds.simpleName}.selection.faa", emit: backbone

  publishDir "${params.output_folder}/backbone", mode: 'copy'


  script:
    """
    $projectDir/gv-phylo/prepare_backbone.py $gvdg_genomes_tsv $seeds ${params.selectors} -i ${params.ingroup} -s ${params.subgroup}
    """
}

process concat_seeds_markers {
  input:
    file marker
    each filtered_seed

  output:
    path "${marker.simpleName}.faa", emit: concat_markers

  publishDir "${params.output_folder}/alignments", mode: 'copy'

  when:
    "${marker.simpleName}" == "${filtered_seed.simpleName}"
    
  script:
    """
    cat ${marker} ${filtered_seed} > ${marker.simpleName}.faa
    """
}

process run_mafft {
  input:
    file faa

  output:
    path "${faa.simpleName}.mafft", emit: alignments

  label "mid_cpu"
  publishDir "${params.output_folder}/alignments", mode: 'copy'
  
  script:
    """
    mafft-einsi --thread ${task.cpus} ${faa} > ${faa.simpleName}.mafft
    """
}

process run_mafft_add {
  input:
    file marker
    each aligned_seed
  
  output:
    path "${marker.simpleName}.add.mafft", emit: alignments

  label "mid_cpu"
  publishDir "${params.output_folder}/alignments", mode: 'copy'

  when:
    "${marker.simpleName}" == "${aligned_seed.simpleName}"
      
  script:
    """
    mafft-einsi --thread ${task.cpus} --add ${marker} $aligned_seed > ${marker.simpleName}.add.mafft
    """
}

process run_divvier {
  input:
    file mafft

  output:
    path "${mafft.simpleName}.divvy", emit: divvied

  publishDir "${params.output_folder}/alignments", mode: 'copy'
  
  script:
    """
    divvier -mincol 10 ${mafft}
    sed '/^>/! s/*/-/g' ${mafft}.divvy.fas > ${mafft.simpleName}.divvy
    """
}

process run_trimal {
  input:
    file divvy

  output:
    path "${divvy.simpleName}.trimal", emit: trimaled

  publishDir "${params.output_folder}/alignments", mode: 'copy'
  
  script:
    """
    remove_all_gap_sequences.py -a ${divvy} -o ${divvy.simpleName}.filtered
    trimal -gt .03 -in ${divvy.simpleName}.filtered -out ${divvy.simpleName}.trimal    
    """
}

process run_iqtree {
  input:
    file trimal

  output:
    path "${trimal.simpleName}.iq.treefile", emit: treefiles
    path "${trimal.simpleName}.iq.*", emit: iqfiles

  label "mid_cpu"
  publishDir "${params.output_folder}/trees", mode: 'copy'
  
  script:
    if( params.treemode == 'iqtree-fast' )
      """
      iqtree -s ${trimal} \
             -m LG+G+I \
             -fast \
             -pre ${trimal.simpleName}.iq \
             -nt ${task.cpus}
      """
    else if( params.treemode == 'iqtree-ufboot' )
      """
      iqtree -s ${trimal} \
             -mset LG -mrate E,I,G,I+G,R \
             -bb 1000 -bnni \
             -pre ${trimal.simpleName}.iq \
             -nt ${task.cpus}
      """
}

process color_tree {
  input:
    file treefile
    file colors

  output:
  path "${treefile.simpleName}.nex", emit: nexus

  publishDir "${params.output_folder}/trees", mode: 'copy'
  
  script:
    """
    color_tree.py \
            -i ${treefile} \
            -c ${colors} \
            -o ${treefile.simpleName}.nex
    """
}
