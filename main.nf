#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//--------------------------------------------------------------------------
// Param Checking
//--------------------------------------------------------------------------

if(!params.fastaSubsetSize) {
  throw new Exception("Missing params.fastaSubsetSize")
}

if(params.inputFilePath) {
  seqs = Channel.fromPath( params.inputFilePath )
           .splitFasta( by:params.fastaSubsetSize, file:true  )
}
else {
  throw new Exception("Missing params.inputFilePath")
}

//--------------------------------------------------------------------------
// Main Workflow
//--------------------------------------------------------------------------

workflow {
  exportpred(seqs)
  gff = exportpred2gff(exportpred.out)
  indexResults(gff.collectFile(), params.outputFileName)
}

process exportpred {
  container = "jbrestel/exportpred:latest"

input:
    path subsetFasta

  output:
    path 'exportpred_results'

  script:
  """
  exportpred --input $subsetFasta --output exportpred_results
  """
}

process exportpred2gff {
  container = 'bioperl/bioperl:stable'

  input:
    path subset

  output:
    path 'exportpred_subset.gff'

  script:
  """
  exportpred2gff.pl --inputFile $subset \
                 --outputFile exportpred_subset.gff
  """

}

process indexResults {
  container = 'biocontainers/tabix:v1.9-11-deb_cv1'

  publishDir params.outputDir, mode: 'copy'

  input:
    path gff
    val outputFileName

  output:
    path '*.gz'
    path '*.gz.tbi'

  script:
  """
  sort -k1,1 -k4,4n $gff > ${outputFileName}
  bgzip ${outputFileName}
  tabix -p gff ${outputFileName}.gz
  """
}
