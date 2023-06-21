/*
 * pipeline input parameters
 */
params.source_folder = "$projectDir/lib/R/"
params.direction = "undirected"
params.BN_filename = "$projectDir/input/BN_digraph_pruned_formatted"
params.genelist_filename = "$projectDir/input/KDAInputFile.txt"
params.outputDir = "results/"
params.layers = 6

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    direction              : ${params.direction}
    layers                 : ${params.layers}
    Bayesian network       : ${params.BN_filename}
    Genelist filename      : ${params.genelist_filename}
    """
    .stripIndent()

process RUN_KDA {
    publishDir params.outputDir, mode:'copy'

    input:
    val source_folder
    val direction
    path BN_filename
    path input_filename
    //val outputDir
    val layers

    script: 
    """
    R-keydriver-analysis.R $source_folder $direction $BN_filename $input_filename $layers 
    """
}

workflow {
    
    runkda_ch = RUN_KDA(params.source_folder,params.direction,params.BN_filename,params.genelist_filename,params.outputDir,params.layers)

}
