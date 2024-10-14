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
    | collectFile(name: params.outputFileName, storeDir: params.outputDir)
}

process exportpred {

input:
    path subsetFasta

  output:
    path 'exportpred_results'

  script:
  """
  exportpred --input $subsetFasta --output exportpred_results
  """
}
