## check if the data is balanced
isbalanced <- function(data)
{
  suppressWarnings(!is.list(replications(~ . , data)))
}

## PER'S FUNCTION MAMANALYSIS.R
runMAM <- function(data, Prod_effects, individual, attributes, adjustedMAM=FALSE, 
                   alpha_conditionalMAM=1){
  if(length(attributes) < 2)
    stop("number of attributes for MAM should be more than 1")
  if(length(Prod_effects) > 1)
    stop("should be one-way product structure")  
  dataMAM <- data[, c(individual, Prod_effects)]
  dataMAM$replication <- rep(0, nrow(data))
  dataMAM[, 1] <- as.factor(dataMAM[, 1])
  dataMAM[, 2] <- as.factor(dataMAM[, 2])
  if(nlevels(dataMAM[, 2]) < 3)
    stop("There MUST be at least 3 products")
  ## create a rep factor
  assprod <- interaction(dataMAM[, 1], dataMAM[, 2])
  t <- table(assprod)
  if(length(unique(t))!=1)
    stop("data is unbalanced")
  for(i in 1:length(names(t)))
    dataMAM$replication[assprod==names(t)[i]] <- 1:unique(t)
  dataMAM <- cbind(dataMAM, data[, attributes])
  return(MAManalysis(dataMAM, adjustedMAM, alpha_conditionalMAM))
}


### function checks  if there are zero cells in a factor term
## UNUSED since version 2.0-7
checkZeroCell <- function(data, factors)
{
  t <- table(data[, match(factors, names(data))])
  if(length(which(t==0))>0)
  {
    message(paste("Some of the combinations of ", paste(factors,collapse=":"), 
                  " has no data, therefore this combination will not be part of the initial model"))
    cat("\n")
    return(TRUE)
  }
  
  return(FALSE)
}

## checks if the number of levels for an interaction term 
## is equal to number of observations
checkNumberInteract <- function(data, factors)
{
  ## returns TRUE if number of levels is equal to nrow of data
  
  nlev <- 1
  for(i in 1:length(factors))
  {    
    if(!is.factor(data[, match(factors[i], names(data))]))
      next()
    nlev <- nlev * nlevels(data[, match(factors[i], names(data))])
  }
  if(nlev >= nrow(data))
  {
    warning.str <- "Number of levels for "
    if(length(factors) > 1)
      warning.str <- c(warning.str," interaction ", sep=" ")
 
    warning.str <- c(warning.str, paste(factors,collapse=":"), 
                     " is more or equal to the number of observations in data", 
                     sep=" ")    
    message(warning.str)
    cat("\n")
    return(TRUE)
  }
  return(FALSE)
}


### Function converts variables to factors
convertToFactors <- function(data, facs)
{
  #convert effects to factors
  for(fac in facs)
    data[,fac] <- as.factor(data[,fac])
  data
}

## create formula with only fixed terms
fixedOrRandFormula <- function(fmodel, isfixed = TRUE)
{
  terms.fm <- attr(terms.formula(fmodel),"term.labels")
  ind.rand.terms <- which(unlist(lapply(terms.fm,
                                        function(x) 
                                          substring.location(x, "|")$first))!=0)
  terms.fm[ind.rand.terms] <- unlist(lapply(terms.fm[ind.rand.terms],
                                            function(x) paste("(",x,")",sep="")))
  fm <- paste(fmodel)
  if(isfixed)
    fm[3] <- paste(terms.fm[-ind.rand.terms], collapse=" + ")
  else
    fm[3] <- paste(terms.fm[ind.rand.terms], collapse=" + ")
  if(fm[3]=="")
    fo <- as.formula(paste(fm[2],fm[1],1, sep=""))
  else
    fo <- as.formula(paste(fm[2],fm[1],fm[3], sep=""))
  return(fo)
}


.createFormulaAnovaLsmeans <- function(mf.final, mf.final.lsm, random, fixed,
                                       mult.scaling, data, oneway_rand = TRUE){
  
  ff <- fixedOrRandFormula(mf.final)
  fr <- fixedOrRandFormula(mf.final.lsm, isfixed = FALSE)
  ## maximal order of interaction
  max.inter <- max(attr(terms(ff), "order")) 
  scaling.private.adds <- NULL
  if(max.inter > 2){
    scaling.private.adds <- paste("scaling.private.add", 1:(max.inter - 2) , sep="")
    for (i in 1:length(scaling.private.adds)){
      assign(scaling.private.adds[i], rep(1, nrow(data)) )
    }
    data[, scaling.private.adds] <- rep(1, nrow(data))     
  }
    
  
  
  if(length(fixed$Product)>1 && mult.scaling){
    prods <- paste("x.scaling.private", fixed$Product, sep="")
    for(namesProd in fixed$Product){
      assign(paste("x.scaling.private", namesProd, sep=""), 
             scale(predict(lm(as.formula(paste(ff[2], ff[1], namesProd)), 
                              data=data)), scale=FALSE))
    }

  }
  
  # create x out of predicted values from lm
  data$x.scaling.private <- rep(NA, nrow(data))
  x.scaling.private <- scale(predict(lm(ff, data=data)), scale=FALSE)
  notNA <- rownames(x.scaling.private)
  data[notNA, "x.scaling.private"] <- x.scaling.private
  
  ## for anova
  fm <- paste(mf.final)
  if(is.list(random))
    rnd.fo <- random$individual
  else
    rnd.fo <- random

  if(length(fixed$Product) > 1 && mult.scaling){
    if(!is.null(scaling.private.adds))
      fm[3] <- paste(fm[3], paste(rnd.fo, prods, 
                                paste(scaling.private.adds, collapse = ":"), 
                                sep=":", 
                                collapse=" + "), sep =" + ")
    else
      fm[3] <- paste(fm[3], paste(rnd.fo, prods, 
                                  sep=":", 
                                  collapse=" + "), sep =" + ")
  }
  else{
    if(oneway_rand)
      rightff <- paste(paste(ff)[3], paste(fr)[3], sep = " + ")
    else
      rightff <- fm[3]
    if(!is.null(scaling.private.adds))
      fm[3] <- paste(rightff, paste(rnd.fo, 
                                "x.scaling.private", 
                                paste(scaling.private.adds, collapse = ":"),
                                sep=":"), sep=" + ")
    else
      fm[3] <- paste(rightff, paste(rnd.fo, 
                                  "x.scaling.private",
                                  sep=":"), sep=" + ")
  }

 
  fo.anova <- as.formula(paste(fm[2], fm[1], fm[3], sep=""))    
  
return(list(fo.anova = fo.anova, data = data))
} 

### Create an lmer model
createLMERmodel <- function(structure, data, response, fixed, random, corr, 
                            MAM = FALSE, mult.scaling = FALSE, 
                            calc_post_hoc = FALSE,  oneway_rand = TRUE)
{ 
  
  #construct formula for lmer model    
  mf.final <- createFormulaAllFixRand(structure, data, response, fixed, random, 
                                      corr)    
  
  
  ## if MAM needs to be constructed
  if(MAM){
    if(length(fixed$Product) > 1){
      data$prod <- interaction(data[, fixed$Product])
      mf.final.lsm <- createFormulaAllFixRand(structure, data, response, 
                                              list(Product="prod", 
                                                   Consumer=fixed$Consumer), 
                                              random, corr)   
    }else{
      mf.final.lsm <- mf.final
    }
    
    
    
    ## create formulas for anova and lsmeans   
    ############################################################################
    fo <- .createFormulaAnovaLsmeans(mf.final, mf.final.lsm, random, fixed,
                                     mult.scaling, data, 
                                     oneway_rand = oneway_rand)
    
    ## for anova
    modelMAM <- lmerTest::lmer(fo$fo.anova, fo$data)

    return(modelMAM = modelMAM)   
  }else{
    model <- lmerTest::lmer(mf.final, data) 
    return(model)
  }
    
}


# check an interaction term for validity
checkComb <- function(data, factors)
{
  ## removed checkZeroCell in 2.0-7 version
  ## since Rune has fixed an issue with rank deficiency
  return(checkNumberInteract(data,factors))# || checkZeroCell(data, factors))
}

.fixedrand <- function(model)
{
  effs <- attr(terms(formula(model)), "term.labels")
  neffs <- length(effs)
  randeffs <- effs[grep(" | ", effs)]
  randeffs <- sapply(randeffs, function(x) substring(x, 5, nchar(x)))
  fixedeffs <- effs[!(effs %in% names(randeffs))]
  return(list(randeffs=randeffs, fixedeffs=fixedeffs))
}

.fillpvalues <- function(x, pvalue)
{
  pvalue[rownames(x$anova.table),x$response] <- x$anova.table[,6]
  pvalue
}

.renameScalingTerm <- function(tableWithScaling, Prod_effects){
  idsScaling <- unlist(lapply(rownames(tableWithScaling), 
                function(x) grepl(":x.scaling.private", x)))
  if(sum(idsScaling) > 1){
    for(i in 1:length(Prod_effects)){
      numscale <- unlist(lapply(rownames(tableWithScaling)[idsScaling], 
                    function(x) grepl(Prod_effects[i], x)))
      rownames(tableWithScaling)[idsScaling][numscale] <- paste("Scaling", 
                                                                Prod_effects[i],
                                                                sep = " ")
    }
  }
  else{
    rownames(tableWithScaling)[unlist(lapply(rownames(tableWithScaling), 
                                             function(x) grepl(":x.scaling.private", x)))] <- 
      "Scaling"
  }
   
  tableWithScaling
}


###############################################################################
# get terms contained  - from lmerTest package
###############################################################################
getIndTermsContained <- function(allterms, ind.hoi)
{
  
  terms.hoi.split <- strsplit(allterms[ind.hoi],":")
  ind.terms.contain <- NULL
  #check which of the terms are contained in the highest order terms
  for(i in (1:length(allterms))[-ind.hoi]) 
  {
    isContained<-FALSE
    for(j in 1:length(terms.hoi.split))
    {
      #if the term is contained in some of the highest order interactions then 
      #we cannot test it for significance
      if(length(which(unlist(strsplit(allterms[i],":")) %in% terms.hoi.split[[j]] == FALSE))==0)
      {
        isContained <- TRUE
        break
      }                
    }
    if(isContained)
      ind.terms.contain <- c(ind.terms.contain,i)
    
  }
  # if there are no terms that are contained in the maximum order effects
  # then compare all the terms between each other for the maximum p value
  if( is.null(ind.terms.contain) )
    return(NULL)
  return(ind.terms.contain)
}

###############################################################################
# get terms contained 
###############################################################################
.getIndTermsContained <- function(allterms, ind.hoi)
{
  
  terms.hoi.split <- strsplit(allterms[ind.hoi],":")
  ind.terms.contain <- NULL
  #check which of the terms are contained in the highest order terms
  for(i in (1:length(allterms))[-ind.hoi]) 
  {
    isContained<-FALSE
    for(j in 1:length(terms.hoi.split))
    {
      #if the term is contained in some of the highest order interactions then 
      #we cannot test it for significance
      if(length(which(unlist(strsplit(allterms[i],":")) %in% terms.hoi.split[[j]] == FALSE))==0)
      {
        isContained <- TRUE
        break
      }                
    }
    if(isContained)
      ind.terms.contain <- c(ind.terms.contain,i)
    
  }
  # if there are no terms that are contained in the maximum order effects
  # then compare all the terms between each other for the maximum p value
  if( is.null(ind.terms.contain) )
    return(NULL)
  return(ind.terms.contain)
}


## get the pure lsmeans for an interaction term
getPureInter <- function(lsm.table, anova.table, eff){

  rows.lsm <- sapply(rownames(lsm.table), 
                     function(x) strsplit(x, " ")[[1]][1]) 
  pure.inter.lsm <- lsm.table[which(rows.lsm %in% eff), ]
  
  contained.effs <- 
    rownames(anova.table)[.getIndTermsContained(rownames(anova.table), 
                                               which(rownames(anova.table) 
                                                     == eff))]
  ## deltas for 3 way interactions
  if( length(unlist(strsplit(eff,":"))) == 3 ){
    ind.inteffs <- grep(":", contained.effs)
    ##plust main effs
    main.effs <- contained.effs[-ind.inteffs]
    ##minus the interactions
    contained.effs <- contained.effs[ind.inteffs]
  }

  p1 <- pure.inter.lsm[ , 1:which(colnames(pure.inter.lsm)=="Estimate"), 
                       drop = FALSE]
  for(ceff in contained.effs){
    p1  <- merge(p1, 
                 lsm.table[rows.lsm == ceff, 
                           c(unlist(strsplit(ceff,":")), "Estimate")], 
                 by = unlist(strsplit(ceff,":")))
    p1[, "Estimate.x"] <- p1[, "Estimate.x"] - p1[, "Estimate.y"]
    colnames(p1)[which(colnames(p1) == "Estimate.x")] <- "Estimate"
    p1 <- p1[ ,- which(colnames(p1) == "Estimate.y")]        
  }
  ## plus the main effs
  if( length(unlist(strsplit(eff,":"))) == 3 ){
    for(ceff in main.effs){
      p1  <- merge(p1, 
                   lsm.table[rows.lsm == ceff, 
                             c(unlist(strsplit(ceff,":")), "Estimate")], 
                   by = unlist(strsplit(ceff,":")))
      p1[, "Estimate.x"] <- p1[, "Estimate.x"] + p1[, "Estimate.y"]
      colnames(p1)[which(colnames(p1) == "Estimate.x")] <- "Estimate"
      p1 <- p1[ ,- which(colnames(p1) == "Estimate.y")]        
    }
  }
  p1
}

## calculate pure diffs
.calcPureDiffs <- function(pureinter){
  puredifs <- matrix(0, ncol=nrow(pureinter), nrow = nrow(pureinter))
  for (i in 1:nrow(pureinter)) for (j in 1:i) 
    puredifs[i,j] <- pureinter[i, "Estimate"] -  pureinter[j, "Estimate"]  
  puredifs
}

## calculate average d prime from the step function
.calcAvDprime <- function(model, anova.table, dlsm.table, lsm.table){
  sigma <- summary(model, "lme4")$sigma
  rows <- sapply(rownames(dlsm.table), 
                 function(x) strsplit(x, " ")[[1]][1]) 
  anova.table$dprimeav <- rep(1, nrow(anova.table))
  num.scale <- which(grepl("Scaling", rownames(anova.table)) == TRUE)
  if(length(num.scale) > 0)
    anova.table.noscaleff <- anova.table[-num.scale,]
  else
    anova.table.noscaleff <- anova.table
    
  
  for(eff in rownames(anova.table.noscaleff)){
    
      
    pureinter <- getPureInter(lsm.table, anova.table.noscaleff, eff)
    puredifs <- .calcPureDiffs(pureinter) 
      
      
    #all.equal(sum(pure.inter.pairs^2)/(12), sum(puredifs^2)/(12), tol = 1e-4)
    dp <- puredifs / sigma      
    av.dp <- sqrt(sum(dp^2)/(nrow(dp)*(nrow(dp)-1)/2))
   
    anova.table[eff, "dprimeav"] <- av.dp 
  }
  anova.table <- anova.table[, 
                             c("Sum Sq", "Mean Sq", "NumDF", "DenDF", "F.value",
                               "dprimeav", "Pr(>F)")]
  anova.table
  
 
}



.stepAllAttrNoMAM <- function(attr, product_structure, error_structure,
                              data, Prod_effects, random,
                              reduce.random = reduce.random, 
                              alpha.random = alpha.random, 
                              alpha.fixed = alpha.fixed, 
                              calc_post_hoc = calc_post_hoc,
                              keep.effs){
  m <- suppressMessages(createLMERmodel(structure =
                                          list(product_structure = 
                                                 product_structure,
                                               error_structure = 
                                                 error_structure), 
                                        data = data, response = attr,
                                        fixed = list(Product = Prod_effects,
                                                     Consumer = NULL),
                                        random = random, corr = FALSE,
                                        calc_post_hoc = calc_post_hoc))
  #m <- elimZeroVar(m)
  suppressMessages(st <- step(m, reduce.fixed = FALSE, 
                             reduce.random = reduce.random, 
                             alpha.random = alpha.random, 
                             alpha.fixed = alpha.fixed, 
                             lsmeans.calc = TRUE,
                             difflsmeans.calc = calc_post_hoc,
                             keep.effs = keep.effs))
  
  if(calc_post_hoc)
    st$anova.table <- .calcAvDprime(st$model, st$anova.table, 
                                    st$diffs.lsmeans.table, st$lsmeans.table)

  st
}


## step function for MAM
.stepAllAttrMAM <- function(attr, product_structure, error_structure,
                            data, Prod_effects, random,
                            reduce.random = reduce.random, 
                            alpha.random = alpha.random, 
                            alpha.fixed = alpha.fixed, 
                            mult.scaling = mult.scaling, 
                            calc_post_hoc = calc_post_hoc,
                            keep.effs, oneway_rand = oneway_rand){
  modelMAM <- suppressMessages(createLMERmodel(structure = 
                                                   list(product_structure =
                                                          product_structure, 
                                                        error_structure =
                                                          error_structure), 
                                                 data = data, response = attr,
                                                 fixed = list(Product = 
                                                                Prod_effects, 
                                                              Consumer=NULL),
                                                 random = random, corr = FALSE,
                                                 MAM = TRUE,
                                                 mult.scaling = mult.scaling, 
                                                 calc_post_hoc = calc_post_hoc,
                                                 oneway_rand = oneway_rand))
  
  modelMAM  <- elimZeroVar(modelMAM)

  
  st <- suppressMessages(step(modelMAM, fixed.calc = FALSE,
                              keep.effs = keep.effs, 
                              reduce.random = reduce.random))
  rand.table <- st$rand.table
  
  if(reduce.random){
    modelMAM <- as(st$model, "merModLmerTest")
    anova.table <- suppressMessages(anova(modelMAM, type = 1))
  }
  else
    anova.table <- suppressMessages(anova(modelMAM, type = 1))
  
  anova.table <- .renameScalingTerm(anova.table, Prod_effects) 
  
  if(calc_post_hoc){    
    model <- .createModelFromMAM(modelMAM)
    rho <- rhoInitLsmeans(model, modelMAM)
    lsm.tab <- calcLSMEANS(model, data, rho, alpha = alpha.fixed, 
                           lsmeansORdiff = TRUE)$summ.data
    dlsm.tab <- calcLSMEANS(model, data, rho, alpha = alpha.fixed, 
                            lsmeansORdiff = FALSE)$summ.data
    anova.table <- .calcAvDprime(modelMAM, anova.table, 
                                 dlsm.tab, lsm.tab)
    return(list(anova.table = anova.table, rand.table = rand.table,
                lsmeans.table = lsm.tab, diffs.lsmeans.table = dlsm.tab)) 
  }
  
  return(list(anova.table = anova.table, rand.table = rand.table)) 
}

.fixedrand <- function(model)
{  
  effs <- attr(terms(formula(model)), "term.labels")
  neffs <- length(effs)
  randeffs <- effs[grep(" | ", effs)]
  randeffs <- sapply(randeffs, function(x) substring(x, 5, nchar(x)))
  fixedeffs <- effs[!(effs %in% names(randeffs))]
  return(list(randeffs=randeffs, fixedeffs=fixedeffs))
}

.getAssessor <- function(random){
  if(is.list(random))
    return(random$individual)
  else
    return(random)
}

.createModelFromMAM <- function(modelMAM){
  fmodelMAM <- formula(modelMAM)
  fm <- paste(fmodelMAM)
  terms.fm <- attr(terms(fmodelMAM), "term.labels")
  
  is.scale <- grepl(":x.scaling", terms.fm)   

  fmodel <- paste(fm[2],fm[1], paste(fm[3], paste(terms.fm[is.scale], 
                                                  collapse = " - "), 
                                     sep = " - "))
  return(updateModel(modelMAM, fmodel, getREML(modelMAM), 
              attr(model.matrix(modelMAM),"contrasts")))

} 

