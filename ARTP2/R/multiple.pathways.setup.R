

multiple.pathways.setup <- function(summary.files, pathway, family, reference, lambda, 
                                    ncases, ncontrols, nsamples, options){
  
  start.time <- date()
  
  # validate the format of main inputs
  validate.summary.input(summary.files, pathway, family, reference, lambda, 
                         ncases, ncontrols, nsamples)
  
  # reformat reference to convert factor to character
  reference <- reformat.reference.path(reference)
  
  # merge and reset options
  options <- options.setup(options, family, lambda, ncases, ncontrols, nsamples)
  
  # load definition of pathway
  pathway <- load.pathway.definition(pathway, options)
  
  # load and check summary statistics
  sum.stat <- load.summary.statistics(summary.files, family, pathway$SNP, options)
  
  # estimate P and SE if they are not provided by users
  sum.stat <- complete.sum.stat(sum.stat, options)
  
  # split super-pathway into independent sub-pathways
  sub.pathway <- create.sub.pathway(pathway)
  
  # split summary statistics into sub-pathways
  sum.stat <- split.sum.stat(sum.stat, sub.pathway)
  
  # load SNPs and their reference and effect alleles in reference genotype
  allele.info <- load.reference.allele(reference, pathway, options)
  
  setup <- create.super.pathway.setup(sum.stat, allele.info, reference, sub.pathway, options)
  
  setup
  
}

