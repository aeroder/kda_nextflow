/*
 * pipeline input parameters
 */
params.source_folder       = "$projectDir/bin/R/"
params.direction           = "undirected"
params.BN_filename         = "$projectDir/input/BN_digraph_pruned_formatted"
params.genelist_filename   = "$projectDir/input/DKAinputfile.txt"
params.outputDir           = "results/"
params.layers              = 6

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    direction    : ${params.transcriptome_file}
    levels       : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()
