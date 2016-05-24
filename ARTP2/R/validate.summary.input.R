
validate.summary.input <- function(summary.files, pathway, family, reference, lambda, 
                                   ncases, ncontrols, nsamples){
  
  if(!is.null(summary.files)){
    validate.summary.files(summary.files)
  }
  
  if(!is.null(pathway)){
    validate.pathway.definition(pathway)
  }
  
  if(!is.null(family)){
    validate.family(family)
  }
  
  if(!is.null(reference)){
    validate.reference(reference)
  }
  
  if(!is.null(summary.files) && !is.null(lambda)){
    validate.lambda.summaryData(summary.files, lambda)
  }
  
  if(!is.null(family) && !is.null(lambda) && !is.null(ncases) && !is.null(ncontrols) && !is.null(nsamples)){
    validate.sample.size(family, lambda, ncases, ncontrols, nsamples)
  }
  
}

