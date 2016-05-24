
###############################################################################
# performs  analysis of sensory data
###############################################################################
sensmixedFun <- function(attributes = NULL, Prod_effects, replication = NULL, 
                         individual, data, product_structure = 3, 
                         error_structure = "No_Rep", 
                         MAM = FALSE, mult.scaling = FALSE,  
                         oneway_rand = TRUE, MAM_PER = FALSE, 
                         adjustedMAM = FALSE, 
                         alpha_conditionalMAM = 1, calc_post_hoc = FALSE, 
                         parallel=FALSE, 
                         reduce.random = TRUE, alpha.random = 0.1, 
                         alpha.fixed = 0.05, interact.symbol = interact.symbol,
                         keep.effs = NULL)
  
{
  ## product_structure=1  (default structure) : Analysis of main fixed effects
  ## product_structure=2 : Main effects AND all 2-factor interactions. 
  ## product_structure=3 : Full factorial model with ALL possible fixed effects
  ## error_structure 
  if(is.null(attributes))
    attributes <- colnames(data)[!sapply(data, 
                                         function(x) "factor" %in% class(x))]
  
  if(class(attributes)=="integer")
    attributes <- colnames(data)[attributes]
  
  if(length(attributes) < 7)
    parallel <- FALSE
  
  
  if(MAM_PER)
  {
      return(runMAM(data, Prod_effects, individual, attributes, adjustedMAM, 
                    alpha_conditionalMAM))   
  }

  
  if(!is.null(replication) && error_structure!="No-Rep")
    random <- list(individual=individual, replication=replication)
  else
    random <- individual
  isRandReduce <- TRUE
  isFixReduce <- FALSE
  isLsmeans <- TRUE
  
  
  
 

  nbvar <- length(attributes)
  
  
  
  ## create the initial model
  suppressMessages(model.init <- createLMERmodel(structure = 
                                  list(product_structure=product_structure, 
                           error_structure = error_structure), 
                           data, attributes[1],
                           fixed = list(Product = Prod_effects, Consumer = NULL),
                           random = random, corr=FALSE, MAM, 
                           mult.scaling = mult.scaling,  
                           calc_post_hoc = calc_post_hoc, 
                           oneway_rand = oneway_rand))
  model <- model.init #if(MAM) model.init$modelMAM else model.init

  
 
  #if(checkCorr(model))
  #  isRandReduce <- FALSE
  
  #number of effects in the model
  fixedrand <- .fixedrand(model)
  # number of fixed effects
  nfixe <- length(fixedrand$fixedeffs)

  # number of random effects
  nrand <- length(fixedrand$randeffs)
  
  n <- nfixe + nrand
  
  # Matrix that will contain the p values
  pvalueF <- matrix(NA, nfixe, nbvar)
  colnames(pvalueF) <- attributes
  pvalueChi <- matrix(NA, nrand, nbvar)
  colnames(pvalueChi) <- attributes
  
  ## TODO: change rownames for pvaues according to elimRand.R
  ##       should work for slopes as well
  rownames(pvalueF) <- fixedrand$fixedeffs
  pvalueF <- .renameScalingTerm(pvalueF, Prod_effects) 
  
  rownames(pvalueChi) <- fixedrand$randeffs

  Fval <- matrix(0, nfixe, nbvar)
  colnames(Fval) <- attributes
  rownames(Fval) <- rownames(pvalueF)
  Chi <- matrix(0, nrand, nbvar)
  colnames(Chi) <- attributes
  rownames(Chi) <- rownames(pvalueChi)
 
  ### using parallel
  ## TODO: FIX the use of parallel calculations
  #if(parallel)
  #  return(NULL)
#  else{    
    if(!MAM){      
#       res <- llply(data[,attributes], .stepAllAttrNoMAM, model, reduce.random,
#                    alpha.random, alpha.fixed, calc_post_hoc, .progress="text")
      if(!parallel)
        res <- llply(attributes, .stepAllAttrNoMAM, 
                     product_structure = product_structure, 
                     error_structure = error_structure,
                     data = data, Prod_effects = Prod_effects, random = random,
                     reduce.random = reduce.random, alpha.random = alpha.random, 
                     alpha.fixed = alpha.fixed,  calc_post_hoc, 
                     keep.effs = keep.effs, .progress="text")
      else
        res <- stepParallelNoMAM(attributes, 
                            product_structure = product_structure, 
                            error_structure = error_structure,
                            data = data, Prod_effects = Prod_effects, 
                            random = random,
                            reduce.random = reduce.random, 
                            alpha.random = alpha.random, 
                            alpha.fixed = alpha.fixed,  calc_post_hoc,
                            keep.effs = keep.effs)
    }
    else{
      if(!parallel)
        res <- llply(attributes, .stepAllAttrMAM, 
                    product_structure = product_structure, 
                    error_structure = error_structure,
                    data = data, Prod_effects = Prod_effects, random = random,
                    reduce.random = reduce.random, alpha.random = alpha.random, 
                    alpha.fixed = alpha.fixed, 
                    mult.scaling = mult.scaling, 
                    calc_post_hoc = calc_post_hoc, 
                    keep.effs = keep.effs, oneway_rand = oneway_rand,
                    .progress="text") 
      else
        res <- stepParallelMAM(attributes,
                   product_structure = product_structure, 
                   error_structure = error_structure,
                   data = data, Prod_effects = Prod_effects, random = random,
                   reduce.random = reduce.random, alpha.random = alpha.random, 
                   alpha.fixed = alpha.fixed, 
                   mult.scaling = mult.scaling, 
                   calc_post_hoc = calc_post_hoc,
                   keep.effs = keep.effs, oneway_rand = oneway_rand)  
    } 
    names(res) <- attributes
 # }
  
  
 
  
  #### fill the results
  finalres <- tryCatch({
    if(calc_post_hoc){
      post_hoc <- vector(mode = "list", length = length(res))      
      names(post_hoc) <- attributes
      ## create d prime output
      dprimeav <- Fval
    }
      
    for(i in 1:length(attributes))
    {    
      res[[i]]$response <- attributes[i]
      
      #fill pvalues for fixed effects
      calc <- res[[i]]$anova.table
      ## if the reduced is lm model
      if("Residuals" %in% rownames(res[[i]]$anova.table))
        pvalueF[rownames(calc)[-nrow(calc)], i] <- calc[-nrow(calc), 5]
      else
        pvalueF[rownames(calc), i] <- calc[, "Pr(>F)"]
      
      #fill pvalues for random effects
      calcrand <- res[[i]]$rand.table    
      pvalueChi[rownames(calcrand),i] <- calcrand[, "p.value"]
      
      # fill F and Chi values
      if("Residuals" %in% rownames(res[[i]]$anova.table))
        Fval[rownames(calc)[-nrow(calc)],i] <- calc[-nrow(calc), 4]
      else
        Fval[rownames(calc),i] <- calc[,"F.value"]
      Chi[rownames(calcrand),i] <- calcrand[,"Chi.sq"] 
      
      ## fill differences of lsmeans
      if(calc_post_hoc){
        post_hoc[[i]] <- res[[i]]$diffs.lsmeans.table
        dprimeav[rownames(calc), i] <- calc[, "dprimeav"]
      }
    }
    
    pvalueF[is.na(pvalueF)] <- 1
    pvalueChi[is.na(pvalueChi)] <- 1  
      

    }, error = function(e) { NULL })
  if(is.null(finalres) && parallel){
    message(" \n WARNING: error in parallel has occurred: cannot call Kenward-Roger 
            in parallel. instead use unparallelized version \n")     
    return(sensmixedFun(attributes = attributes, Prod_effects = Prod_effects, 
                        replication = replication, individual = individual, 
                        data = data, product_structure = product_structure, 
                        error_structure = error_structure, MAM = MAM, 
                        mult.scaling = mult.scaling,  
                        oneway_rand = oneway_rand, 
                        MAM_PER = MAM_PER, adjustedMAM = adjustedMAM, 
                        alpha_conditionalMAM = alpha_conditionalMAM, 
                        calc_post_hoc = calc_post_hoc, parallel = FALSE, 
                        reduce.random = reduce.random, 
                        alpha.random = alpha.random, alpha.fixed = alpha.fixed, 
                        interact.symbol = interact.symbol,
                        keep.effs = keep.effs))
  }
  ## change the output
  ## output for the random effects
  #tr_rand <- .changeOutput(Chi, pvalueChi, TRUE)
    
  
  if(MAM){   
    ind.scaling <- grepl("Scaling", rownames(pvalueF))
    pvalueScaling <- pvalueF[ind.scaling, , drop=FALSE]
    pvalueF <- pvalueF[!ind.scaling, , drop=FALSE]
    FScaling <- Fval[ind.scaling, , drop=FALSE]
    Fval <- Fval[!ind.scaling, , drop=FALSE]  
    
    #tr_fixed <- .changeOutput(Fval, pvalueF, FALSE)
    #tr_scaling <- .changeOutput(FScaling, pvalueScaling, FALSE)    
    if(calc_post_hoc){
      dprimeav <- dprimeav[!ind.scaling, , drop=FALSE]
      return(list(fixed = list(Fval = Fval, dprimeav = dprimeav, 
                               pvalueF = pvalueF), 
                random = list(Chi = Chi, pvalueChi = pvalueChi), 
                scaling = list(FScaling = FScaling, 
                              pvalueScaling = pvalueScaling), 
                post_hoc = post_hoc, step_res = res))
    }
    else
      return(list(fixed = list(Fval = Fval, pvalueF = pvalueF), 
                  random = list(Chi = Chi, pvalueChi = pvalueChi), 
                  scaling = list(FScaling = FScaling, 
                                 pvalueScaling = pvalueScaling),
                  step_res = res))
  }
  #tr_fixed <- .changeOutput(Fval, pvalueF, FALSE)
  
  if(calc_post_hoc)    
      return(list(fixed = list(Fval = Fval, dprimeav = dprimeav, 
                               pvalueF = pvalueF), 
                  random = list(Chi = Chi, pvalueChi = pvalueChi), 
                  post_hoc = post_hoc, step_res = res))  
  return(list(fixed = list(Fval = Fval, pvalueF = pvalueF), 
              random = list(Chi = Chi, pvalueChi = pvalueChi), step_res = res))  
}
