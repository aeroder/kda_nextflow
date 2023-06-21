log.info """\
    K D A - N F   P I P E L I N E
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
    //val source_folder
    val direction
    path BN_filename
    path input_filename
    //val outputDir
    val layers

    output:
    path '*keydriver.txt'
    path '*combined.txt'
    path '*parameters.txt'

    script: 
    """
    R-keydriver-analysis.R $direction $BN_filename $input_filename $layers $projectDir/bin/functions.R  
    """
}

workflow {

    runkda_ch = RUN_KDA(params.direction,params.BN_filename,params.genelist_filename,params.layers)

}
