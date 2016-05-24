

multiple.comparison.summaryData <- function(summary.files, pathway, family, reference, lambda, 
                                            ncases = list(), ncontrols = list(), nsamples = list(), 
                                            options = NULL){
  
  options(warn = 1)
  
  validate.summary.input(summary.files, NULL, family, reference, lambda, ncases, ncontrols, nsamples)
  
  reference <- reformat.reference.path(reference)
  
  options <- options.setup(options, family, lambda, ncases, ncontrols, nsamples)
  
  pathway <- load.pathway.set(pathway, options)
  
  super.pathway <- create.super.pathway(pathway)
  
  setup <- multiple.pathways.setup(summary.files, super.pathway, family, reference, lambda, 
                                   ncases, ncontrols, nsamples, options)
  
  setup <- recreate.pathway(setup, pathway)
  
  save(setup, file = options$path.setup)
  msg <- paste0("setup file has been saved at ", options$path.setup)
  message(msg)
  
  ret <- generate.pathway.pvalue.stat(setup)
  
  options(warn = 0)
  
  ret
  
}





