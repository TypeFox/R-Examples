
sARTP <- function(summary.files, pathway, family, reference, lambda, 
                                ncases = list(), ncontrols = list(), nsamples = list(), 
                                options = NULL){
  
  options(warn = 1)
  
  setup <- summaryData.setup(summary.files, pathway, family, reference, lambda, 
                             ncases, ncontrols, nsamples, options)
  
  if(setup$options$only.setup){
    return(setup)
  }
  
  test <- norm.stat.test(setup)
  
  options(warn = 0)
  
  list(pathway.pvalue = test$pathway.pvalue, gene.pvalue = test$gene.pvalue, 
       model = test$model, most.sig.genes = test$most.sig.genes, arr.rank = test$arr.rank, 
       accurate = test$accurate, test.timing = test$test.timing, 
       pathway = setup$pathway, deleted.snps = setup$deleted.snps, 
       deleted.genes = setup$deleted.genes, 
       options = setup$options, setup.timing = setup$setup.timing, 
       setup = setup)
  
}



