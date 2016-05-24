
rARTP <- function(formula, data, pathway, family, geno.files = NULL, lambda = 1.0, subset = NULL, options = NULL){
  
  options(warn = 1)
  
  it <- input.type(geno.files)
  
  if(it == 'data.frame'){
    setup <- rawData.dataframe.setup(formula, data, pathway, family, lambda, subset, options)
  }else{
    if(it == 'geno.files'){
      setup <- rawData.genofiles.setup(formula, data, pathway, family, geno.files, lambda, subset, options)
    }else{
      setup <- rawData.plinkfiles.setup(formula, data, pathway, family, geno.files, lambda, subset, options)
    }
  }
  
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



