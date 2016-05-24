
# lambda is the second round adjusted inflation factor
# lambda is different with the ones in ARTP and sARTP, in which lambda is the study-specific first round inflation factor
warm.start <- function(setup, nperm = NULL, lambda = 1.0, nthread = NULL){
  
  options(warn = 1)
  
  setup <- validate.setup(setup)
  
  setup <- update.setup(setup, nperm, lambda, nthread)
  
  test <- norm.stat.test(setup)
  
  options(warn = 0)
  
  list(pathway.pvalue = test$pathway.pvalue, gene.pvalue = test$gene.pvalue, 
       model = test$model, most.sig.genes = test$most.sig.genes, arr.rank = test$arr.rank, 
       accurate = test$accurate, test.timing = test$test.timing, 
       pathway = setup$pathway, meta.stat = setup$meta.stat, 
       deleted.snps = setup$deleted.snps, deleted.genes = setup$deleted.genes, 
       options = setup$options, setup.timing = setup$setup.timing)
  
}
