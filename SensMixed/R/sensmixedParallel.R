parallel <- FALSE

if(parallel){
  respar <- tryCatch({
    funcNOMAM <- local({
      refit
      step
      model
      data
      reduce.random
      alpha.random
      alpha.fixed
      calc_post_hoc
      summary
      
      function(new.resp.private.sensmixed)
      {
        assign("new.resp.private.sensmixed", new.resp.private.sensmixed, 
               envir=environment(formula(model)))
        suppressMessages(m <- refit(object=model, 
                                    newresp=new.resp.private.sensmixed, 
                                    rename.response = TRUE))
        suppressMessages(s <- step(m, reduce.fixed = FALSE, 
                                   reduce.random = reduce.random, 
                                   alpha.random = alpha.random, 
                                   alpha.fixed = alpha.fixed, 
                                   lsmeans.calc=TRUE,
                                   difflsmeans.calc = TRUE))        
        
        if(calc_post_hoc){
 
        }       
        
        s        
      }})
    funcMAM <- local({

      
      function(attr)
      {
        model.init <- suppressMessages(createLMERmodel(
          structure = list(product_structure=product_structure,
                           error_structure=error_structure), 
          data = data, response = attr, fixed = list(Product = Prod_effects, 
                                                     Consumer = NULL),
          random = random, corr = FALSE, MAM = TRUE, 
          mult.scaling = FALSE, 
          calc_post_hoc = calc_post_hoc))
        
        model.an <- model.init$model.anova
        model.lsm <- model.init$model.lsmeans
        #return(model.an)
        
        st <- suppressMessages(step(as(model.an,"merModLmerTest"), 
                                    fixed.calc=FALSE))
        rand.table <- st$rand.table        
        
        if(reduce.random){
          anova.table <- suppressMessages(anova(as(st$model, "merModLmerTest"), 
                                                type = 1))
          if(length(which(anova.table[, "Pr(>F)"] == "NaN") > 0))
            anova.table <- suppressMessages(anova(as(st$model, "merModLmerTest"), 
                                                  type = 1, ddf="Kenward-Roger")) 
        }
        else{
          anova.table <- suppressMessages(anova(model.an, type = 1))
          if(length(which(anova.table[, "Pr(>F)"] == "NaN") > 0))
            anova.table <- suppressMessages(anova(as(model.an, "merModLmerTest"), 
                                                  type = 1, ddf="Kenward-Roger")) 
        }
        
        anova.table <- .renameScalingTerm(anova.table, Prod_effects) 
        
        
        if(length(Prod_effects) > 1)
          lsmeans.table <- lsmeans::.old.lsmeans(model.lsm, pairwise ~ prod)
        else 
          lsmeans.table <-  eval(substitute(lsmeans::.old.lsmeans(object = model.lsm, 
                                                                  pairwise ~ prod), 
                                            list(prod=as.name(Prod_effects))))  
        return(list(anova.table = anova.table, rand.table = rand.table,
                    lsmeans.table = lsmeans.table)) 
       
      }})
    ncpus <- detectCores()
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
                                          "bdiag"
    ),
    #"model", "model.lsm", "data" ),
    envir=environment())
    if(RNGkind()[1L] == "L'Ecuyer-CMRG")
      parallel::clusterSetRNGStream(cl)
    if(!MAM)
      res <- parallel::parLapply(cl, attributes, .stepAllAttrNoMAM, 
                                 product_structure = product_structure, 
                                 error_structure = error_structure,
                                 data = data, Prod_effects = Prod_effects, 
                                 random = random,
                                 reduce.random = reduce.random, 
                                 alpha.random = alpha.random, 
                                 alpha.fixed = alpha.fixed,  calc_post_hoc)
    else
      res <- parallel::parLapply(cl, attributes, funcMAM)
    parallel::stopCluster(cl) 
    res
  },  error = function(e) { NULL })
  if(is.null(respar)){
    message(" \n WARNING: error in parallel has occurred, an unparallized version is used instead \n")    
    return(sensmixedFun(attributes = attributes, Prod_effects = Prod_effects, 
                        replication = replication, individual = individual, 
                        data = data, product_structure = product_structure, 
                        error_structure = error_structure, MAM = MAM, 
                        MAM_PER = MAM_PER, adjustedMAM = adjustedMAM, 
                        alpha_conditionalMAM = alpha_conditionalMAM, 
                        calc_post_hoc = calc_post_hoc, parallel = FALSE, 
                        reduce.random = reduce.random, 
                        alpha.random = alpha.random, alpha.fixed = alpha.fixed, 
                        interact.symbol = interact.symbol))
  }
}