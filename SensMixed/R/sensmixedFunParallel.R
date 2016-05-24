stepParallelNoMAM <- function(attributes, 
                         product_structure = product_structure, 
                         error_structure = error_structure,
                         data = data, Prod_effects = Prod_effects, 
                         random = random,
                         reduce.random = reduce.random, 
                         alpha.random = alpha.random, 
                         alpha.fixed = alpha.fixed,  calc_post_hoc,
                         keep.effs){

  ncpus <- getOption("mc.cores", 2L)# 3#parallel::detectCores()
  cl <- parallel::makePSOCKcluster(ncpus) 
  #cl <- parallel::makeCluster(ncpus) 
  parallel::clusterExport(cl, varlist=c("refit", "data", "createLMERmodel",
                                        "createFormulaAllFixRand",
                                        "checkComb",
                                        "checkNumberInteract",
                                        "checkZeroCell",
                                        "nrandTerms",
                                       # "fixedFormula",
                                        "substring.location", "step",
                                        ".renameScalingTerm",
                                        "isLMM",
                                        "fixef", "anova", 
                                        ".calcAvDprime",
                                        "getPureInter",
                                        ".calcPureDiffs",
                                        ".getIndTermsContained",
                                        "bdiag", "elimZeroVar", "VarCorr"
  ),
  #"model", "model.lsm", "data" ),
  envir=environment())
  if(RNGkind()[1L] == "L'Ecuyer-CMRG")
    parallel::clusterSetRNGStream(cl)
  tryCatch(res <- parallel::parLapply(cl, attributes, .stepAllAttrNoMAM, 
                                      product_structure = product_structure, 
                                      error_structure = error_structure,
                                      data = data, Prod_effects = Prod_effects, 
                                      random = random,
                                      reduce.random = reduce.random, 
                                      alpha.random = alpha.random, 
                                      alpha.fixed = alpha.fixed,  
                                      calc_post_hoc, keep.effs = keep.effs), 
           error = function(e) { NULL }, finally = parallel::stopCluster(cl))
  
return(res)
}


stepParallelMAM <- function(attributes, 
                            product_structure = product_structure, 
                            error_structure = error_structure,
                            data = data, Prod_effects = Prod_effects, random = random,
                            reduce.random = reduce.random, alpha.random = alpha.random, 
                            alpha.fixed = alpha.fixed, 
                            mult.scaling = mult.scaling, 
                            calc_post_hoc = calc_post_hoc,
                            keep.effs = keep.effs, oneway_rand = oneway_rand){
  
  ncpus <- getOption("mc.cores", 2L) #parallel::detectCores()
  cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))             
  parallel::clusterExport(cl, varlist=c("refit", "data", "createLMERmodel",
                                        "createFormulaAllFixRand",
                                        "checkComb",
                                        "checkNumberInteract",
                                        "checkZeroCell",
                                        "nrandTerms",
                                        #"fixedFormula",
                                        "substring.location", "step",
                                        ".renameScalingTerm",
                                        "isLMM",
                                        "fixef", "anova", 
                                        ".calcAvDprime",
                                        "getPureInter",
                                        ".calcPureDiffs",
                                        ".getIndTermsContained",
                                        "bdiag", "elimZeroVar", "VarCorr",
                                        "sigma", "getME", "isGLMM"
  ),
  #"model", "model.lsm", "data" ),
  envir=environment())
  if(RNGkind()[1L] == "L'Ecuyer-CMRG")
    parallel::clusterSetRNGStream(cl)
  tryCatch(res <- parallel::parLapply(cl, attributes, .stepAllAttrMAM, 
                                      product_structure = product_structure, 
                                      error_structure = error_structure,
                                      data = data, Prod_effects = Prod_effects, 
                                      random = random,
                                      reduce.random = reduce.random, 
                                      alpha.random = alpha.random, 
                                      alpha.fixed = alpha.fixed, 
                                      mult.scaling = mult.scaling, 
                                      calc_post_hoc = calc_post_hoc,
                                      keep.effs = keep.effs, oneway_rand = 
                                        oneway_rand), 
           error = function(e) { NULL }, finally = parallel::stopCluster(cl))
  
  return(res)
}