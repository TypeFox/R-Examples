
##########################################################################
getY <- function(model)
{
    return(getME(model, "y"))
}
getX <- function(model)
{
    return(getME(model, "X"))
}
getZt <- function(model)
{
    Ztlist <- getME(model, "Ztlist")#model@Zt
    return(do.call(rBind,Ztlist))
}
getST <- function(model)
{
    return(getME(model, "ST"))
}


##########################################################################
# Create rho environmental variable of mixed model ####################### 
##########################################################################
rhoInitJSS <- function(model)
{
  # creating rho
  rho <- new.env(parent = emptyenv()) # create an empty environment
  rho$model <- model
  rho$y <- getY(model) #model@y                   # store arguments and derived values
  rho$X <- getX(model) #model@X
  chol(rho$XtX <- crossprod(rho$X))       # check for full column rank
  
  rho$REML <-  getREML(model)#model@dims['REML']
  if(class(model) == "lmerMod") 
    rho$s <- summary(model)#, ddf="lme4")
  
  rho$fixEffs <- fixef(model)
  rho$sigma <- sigma(model)
  rho$vlist <- sapply(model@cnms, length)
  
  
  ## get the optima
  ## To UNCOMMENT
  #pp <- model@pp$copy()  
  #opt <- Cv_to_Vv(pp$theta, n = vlist, s = rho$sigma)
  #rho$opt <- opt
  rho$thopt <- getME(model, "theta")
  rho$param <- as.data.frame(VarCorr(model))[, "sdcor"]
#   rho$vars <- Cv_to_Vv(rho$thopt, n = rho$vlist, 
#                        s = rho$sigma)
  return(rho)  
}

       

##############################################################################################
# function to calculate summary of F test with Satterthwaite's approximation of denominator df 
##############################################################################################
calcSatterthJSS  <-  function(Lc, rho)
{
  # F statistics for tested term
  vcov.final <- as.matrix(vcov(rho$model))
  if(is.vector(Lc))
    C.theta.optim <- as.matrix(t(Lc) %*% vcov.final %*% Lc)    
  else
    C.theta.optim <- as.matrix(Lc %*% vcov.final %*% t(Lc))    
  
  invC.theta <- tryCatch({solve(C.theta.optim)}, error = function(e) { NULL })
  if(is.null(invC.theta))
    return(list(denom = 0, Fstat = NA, pvalue = NA, ndf=NA, ss = NA, ms = NA))
  
  q <- qr(C.theta.optim)$rank
  F.stat <- (t(Lc %*% rho$fixEffs) %*% invC.theta %*% (Lc %*% rho$fixEffs))/q
  
  
  #df for F statistics for tested term
  svdec <- eigen(C.theta.optim) 
  
  
  PL <- t(svdec$vectors) %*% Lc
  
  ## based on theta parameters and sigma
  ## also correct
  vss2 <- vcovJSStheta2.temp(rho$model)
  # based on var cor parameters
  # vss <- vcovJSStheta2.var(rho$model)
  theopt <- c(rho$thopt, rho$sigma)
  g <- mygrad(function(x)  vss2(x), theopt)

  if(class(g) == "numeric")
    mat.grad <- llply(1:length(theopt), function(x) matrix(g[x], 
                                                           ncol = ncol(vcov.final), 
                                                           nrow = nrow(vcov.final)))
  else
    mat.grad <- llply(1:length(theopt), function(x) matrix(g[, x], 
                                                         ncol = ncol(vcov.final), 
                                                        nrow = nrow(vcov.final)))
  
  nu.m.fun <- function(m){    
    den.nu <- unlist(llply(1:length(mat.grad), function(x) 
      as.matrix(t(PL[m,]) %*% mat.grad[[x]] %*% PL[m,])))   
    2*(svdec$values[m])^2/(t(den.nu) %*% rho$A %*% den.nu)
  }


  nu.m <- unlist(llply(1:length(svdec$values), .fun = nu.m.fun))

  nu.m[which(abs(2 - nu.m) < 1e-5)] <- 2.00001


  E <- sum( (nu.m/(nu.m-2)) * as.numeric(nu.m>2))
  nu.F <- 2 * E * as.numeric(E > q) / (E - q)
  
  pvalueF <- 1 - pf(F.stat,qr(Lc)$rank, nu.F)
  
  # calculate ss and ms
  if(is.na(F.stat))
    ms <- ss <- NA
  else{
    ms <- F.stat * rho$sigma^2
    ss <- ms * q
  }
  
  ## calculate ss from camp method proc glm
  ## ss <- getSS(Lc, rho$fixEffs ,ginv(rho$XtX)) 
  return( list(ss = ss, ms = ms, denom = nu.F, Fstat = F.stat, 
               pvalue = pvalueF, ndf=q))  
}


## sum of squares based on the proc glm SAS
getSS <- function(L, coef, XtX.) {
  L.beta <- L %*% coef
  if(is.vector(L))
    var.L.beta <- t(L) %*% XtX. %*% L
  if(is.matrix(L))
    var.L.beta <- L %*% XtX. %*% t(L)
  ss <- c(t(L.beta) %*% ginv(var.L.beta) %*% L.beta)
  ss
}

       
###########################################################################
# function to calculate F stat and pvalues for a given term
###########################################################################
calcFpvalueSS <- function(term, Lc, fullCoefs, X.design, model, rho, ddf, 
                           type)
{
  
  if(is.null(Lc))
    return(NULL) 
  
  ## BUG: check vases example from Per
  #calculate ss
  #ss = 1##getSS(Lc, fullCoefs, ginv(crossprod(X.design)))
  {
       
    # for running rune's vcov function
    if(is.vector(Lc))
    {      
      Lc <- Lc[rho$nums.Coefs]
    }  
    else
    {      
      Lc <- Lc[ , rho$nums.Coefs]
    }
  }   
  
   
  if( ddf=="Kenward-Roger" )
  {
    if (!requireNamespace("pbkrtest", quietly = TRUE)) 
      stop("pbkrtest package required for Kenward-Roger's approximations")
    if(is.vector(Lc))
      res.KR <- pbkrtest::KRmodcomp( model, t(as.matrix(Lc)) )
    else
      res.KR <- pbkrtest::KRmodcomp( model, Lc )
    
    ## calculate ms and ss
    ms <- res.KR$test[1,"stat"] * rho$sigma^2
    ss <- ms * res.KR$test[1,"ndf"]
    ## ss <- getSS(Lc, rho$fixEffs ,ginv(rho$XtX)) 
 
    return( list(denom = res.KR$test[1,"ddf"], Fstat = res.KR$test[1,"stat"], 
                 pvalue =  res.KR$test[1,"p.value"], ndf = res.KR$test[1,"ndf"], 
                 ss = ss , ms = ms))
  }
  else
  {
    ## apply satterthwaite's approximation of ddf
    return( c(calcSatterthJSS(Lc, rho)))    
  }
}

###########################################################################
# function to calculate F stat and pvalues for a given term. MAIN
###########################################################################
calcFpvalueMAIN <- function(term, L, X.design, fullCoefs, model, rho, ddf, 
                            type)
{
                         
    if( type == 3 )
    {
      
      Lc <- makeContrastType3SAS(model, term, L)    
      #non identifiable because of rank deficiency
      if(!length(Lc))
        result.fstat <- list(denom=0, Fstat=NA, pvalue=NA, ndf=NA, ss = NA, 
                             ms = NA)
      else
        result.fstat <- calcFpvalueSS(term, Lc, fullCoefs, X.design, model, rho, 
                                      ddf, type)           
    }
    
    if( type == 1 )
    {
      
      find.term <- which(colnames(X.design) == term)
      Lc <- L[find.term[which(find.term %in% rho$nums.Coefs)],]            
      result.fstat <- calcFpvalueSS(term, Lc, fullCoefs, X.design, model, rho, 
                                    ddf, type)      
    } 
    if(type == 2){
      Lc <- makeContrastType2(model, term, L, X.design, rho, fullCoefs) 
      result.fstat <- calcFpvalueSS(term, Lc, fullCoefs, X.design, model, rho, 
                                    ddf, type)
    }
      
   

   c(result.fstat,list(name=term)) 
}


###############################################################################
# function to calculate T test JSS
###############################################################################
calculateTtestJSS <- function(rho, Lc, nrow.res, ddf = "Satterthwaite")
{

  resultTtest <- matrix(0, nrow = nrow.res, ncol = 4)
  colnames(resultTtest) <- c("df", "t value", "p-value", "sqrt.varcor")
  
  if(ddf == "Kenward-Roger"){
    if (!requireNamespace("pbkrtest", quitly = TRUE)) 
      stop("pbkrtest package required for Kenward-Roger's approximations")
    Va <- pbkrtest::vcovAdj(rho$model)
  }
  else{
    # based on theta parameters
    vss <- vcovJSStheta2(rho$model)
    # based on variance parameters
    #vss <- vcovJSStheta2.var(rho$model)
  }
    
  for(i in 1:nrow.res)
  {
    
    if(ddf == "Kenward-Roger"){
      L <- Lc[,i]
      .ddf <- pbkrtest::get_ddf_Lb(rho$model, L)      
      b.hat <- rho$fixEffs
      Lb.hat <- sum(L * b.hat)
      Va.Lb.hat <- t(L) %*% Va %*% L
      t.stat <- as.numeric(Lb.hat / sqrt(Va.Lb.hat))
      p.value <- 2 * pt(abs(t.stat), df = .ddf, lower.tail = FALSE)
      resultTtest[i,1] <- .ddf
      resultTtest[i,2] <- t.stat
      resultTtest[i,3] <- p.value
      resultTtest[i,4] <- as.numeric(sqrt(Va.Lb.hat))
    }
    else{
      ## based on theta parameters
      g <- mygrad(function(x)  vss(t(Lc[,i]), x), c(rho$thopt, rho$sigma))
      ## based on var cor parameters
      #g <- grad(function(x)  vss(t(Lc[,i]), x), rho$vars)
      
      #denominator df
      denom <- t(g) %*% rho$A %*% g
      ## for the theta and sigma parameters
      varcor <- vss(t(Lc[,i]), c(rho$thopt, rho$sigma))
      ## for var cor parameters
      #varcor <- vss(t(Lc[,i]), rho$vars)
      #df
      resultTtest[i,1] <- 2*(varcor)^2/denom
      #statistics
      resultTtest[i,2] <- (Lc[,i] %*%rho$fixEffs)/sqrt(varcor) 
      resultTtest[i,3] <- 2*(1 - pt(abs(resultTtest[i,2]), df = resultTtest[i,1]))
      resultTtest[i,4] <- sqrt(varcor) 
    }
   
  }
  
  return(resultTtest)
}

###############################################################################
# construct design matrix for F test 
###############################################################################
createDesignMat <- function(model, data)
{
model.term <- terms(model)
fixed.term <- attr(model.term,"term.labels") 
X.design <- names.design <-  names.design.withLevels <- NULL

for(i in 1:length(fixed.term))
{

   formula.term <- as.formula(paste("~", fixed.term[i], "- 1"))
   X.design <- cbind(X.design, model.matrix(formula.term, data))
   names.design <- c(names.design, 
                     rep(fixed.term[i],ncol(model.matrix(formula.term, data))))
     
}

if(attr(model.term, "intercept") != 0){
  names.design.withLevels <- c("(Intercept)", colnames(X.design))
  X.design <- cbind(rep(1,dim(X.design)[1]),X.design)
  names.design <- c("(Intercept)", names.design)
}
else
  names.design.withLevels <- colnames(X.design)

colnames(X.design) <- names.design
return(list(X.design=X.design, names.design.withLevels=names.design.withLevels))
}


###############################################################################
# initialize anova table for F test 
###############################################################################
initAnovaTable <- function(model, test.terms, isFixReduce)
{
  anova.table <- matrix(NA, nrow=length(test.terms), ncol=6)
  rownames(anova.table) <- test.terms
  colnames(anova.table) <- c("Sum Sq", "Mean Sq", "NumDF", "DenDF","F.value", 
                             "Pr(>F)")
  anm <- anova(model, ddf="lme4")
  colnames(anm) <- c("NumDF", "Sum Sq", "Mean Sq", "F.value")
  
 
  anova.table[rownames(anm), c("NumDF", "Sum Sq", "Mean Sq")] <- 
    as.matrix(anm[, c("NumDF", "Sum Sq", "Mean Sq")])
  
  if(isFixReduce)
  {
    if(nrow(anova.table)==1)
    {
      anova.table <- c(anova.table[,1:5], 0, anova.table[,6])
      anova.table <- matrix(anova.table, nrow=1, ncol=length(anova.table))
      colnames(anova.table) <- c("Sum Sq", "Mean Sq", "NumDF", "DenDF", 
                                 "F.value", "elim.num","Pr(>F)")
      rownames(anova.table) <- rownames(anm)
      return(anova.table)
    }
    elim.num <- rep(0, nrow(anova.table))
    anova.table <- cbind(anova.table[,1:5], elim.num, anova.table[,6])
    colnames(anova.table)[7] <- "Pr(>F)"
  }
  return(anova.table)    
}


###############################################################################
# get terms contained 
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
      if(length(which(unlist(strsplit(allterms[i],":")) %in% 
                        terms.hoi.split[[j]] == FALSE))==0)
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

orderterms <- function(anova.table)
{
   return(unlist(lapply(rownames(anova.table), function(x) 
     length(unlist(strsplit(x,":"))))))
}
###############################################################################
# get terms to compare in anova.table
###############################################################################
#getTermsToCompare <- function(model)
getTermsToCompare <- function(anova.table, keep.effs = NULL)
{
  
  #order.terms <- attr(terms(model),"order")
  #allterms <- attr(terms(model),"term.labels")
  anova.table.upd <- anova.table[complete.cases(anova.table), , drop=FALSE]
  order.terms <- orderterms(anova.table.upd)
  allterms <- rownames(anova.table.upd )
  ind.hoi <- which(order.terms == max(order.terms))
  ind.terms.contain <- getIndTermsContained(allterms, ind.hoi)
  
  #get the rest of the terms to compare
  allterms.rest <- allterms[-c(ind.terms.contain, ind.hoi)]
  if( length(allterms.rest)==0 )
    terms.compare <- allterms[ind.hoi]
  else
  {
    #get highest order terms in the remaining ones
    order.rest <- unlist(lapply(allterms.rest, function(x) 
      length(unlist(strsplit(x,":")))))
    ind.hoi.rest <- which(order.rest == max(order.rest))
    gtc <- getIndTermsContained(allterms.rest, ind.hoi.rest)
    if( !is.null(gtc) )
      terms.compare <- c(allterms[ind.hoi], 
                         allterms.rest[-getIndTermsContained(allterms.rest, 
                                                             ind.hoi.rest)])
    else
      terms.compare <- c(allterms[ind.hoi], allterms.rest)
  }
  
  return(setdiff(terms.compare, keep.effs))
}

###############################################################################
# find NS effect from the model (starting from highest order interactions)
###############################################################################
getNSFixedTerm <- function(model, anova.table, data, alpha, keep.effs = NULL)
{
  
  pv.max <- 0
  
   if(length(which(anova.table[,"elim.num"]==0))==1)
    terms.compare <- rownames(anova.table)[anova.table[,"elim.num"]==0]
  else
    terms.compare <- getTermsToCompare(anova.table[anova.table[,"elim.num"]==0,], 
                                       keep.effs)
  
  for(tcmp in terms.compare)
  {
    if((!tcmp %in% rownames(anova.table)) || is.na(anova.table[tcmp,"Pr(>F)"]))
      next
    ind <- which(rownames(anova.table)==tcmp)
    if(anova.table[ind, which(colnames(anova.table)=="Pr(>F)")]>=pv.max)
    {
      ns.term <- tcmp
      pv.max <- anova.table[ind,which(colnames(anova.table)=="Pr(>F)")]
    }
  }  
  if(pv.max >= alpha)
    return(ns.term)  
  else
    return(NULL)  
}
  


###############################################################################
# eliminate NS fixed effect from the model
############################################################################### 
elimNSFixedTerm <- function(model, anova.table, data, alpha, elim.num, 
                            l.lmerTest.private.contrast, keep.effs = NULL)
{
  ns.term <- getNSFixedTerm(model, anova.table, data, alpha, 
                            keep.effs = keep.effs)
  if( is.null(ns.term) )
    return(NULL)
  anova.table[ns.term, "elim.num"] <- elim.num
  fm <- formula(model)
  fm[3] <- paste(fm[3], "-", ns.term)
  
  mf.final <- as.formula(paste(fm[2], fm[1], fm[3], sep=""))
  
  model <- updateModel(model, mf.final, getME(model, "is_REML"), 
                       l.lmerTest.private.contrast) 
 
  return( list(model=model, anova.table=anova.table) )
}
  
 
  
#################################################################
# find which effect contains effect term
#################################################################
relatives <- function(classes.term, term, names, factors)
{
  # checks if the terms have the same number of covariates (if any)
  checkCovContain <- function(term1, term2)
  {        
    num.numeric <- which(classes.term=="numeric")
    num.numeric.term1 <- which((num.numeric %in% which(factors[,term1]!=0))==TRUE)
    num.numeric.term2 <- which((num.numeric %in% which(factors[,term2]!=0))==TRUE)
    if((length(num.numeric.term1)>0 && length(num.numeric.term2)>0)||
         (length(num.numeric.term1)==0 && length(num.numeric.term2)==0))
       return(all(num.numeric.term2 == num.numeric.term1))
    else
       return(FALSE)
  }
  is.relative <- function(term1, term2) 
  {
    return(all(!(factors[,term1]&(!factors[,term2]))) && checkCovContain(term1,term2))
  }
  if(length(names) == 1) return(NULL)
  	 which.term <- which(term==names)
	  (1:length(names))[-which.term][sapply(names[-which.term], 
		  			function(term2) is.relative(term, term2))]
}

############################################################################       
# caclulate the General contrast matrix for the hypothesis (as in SAS)
############################################################################
calcGeneralSetForHypothesis <- function(X.design, rho)
{
    
  xtx <- t(X.design) %*% X.design
  
  g2 <- matrix(0,ncol=ncol(xtx), nrow=nrow(xtx))
  
  
  inds <- rho$nums.Coefs
  g2[inds,inds] <- solve(xtx[inds,inds])
  g2[abs(g2)<1e-10] <- 0
  
  
  #general set of estimable function
  L <- g2 %*% xtx
  L[abs(L)<1e-6] <- 0
  return(L)
}
     
       
############################################################################       
# type 3 hypothesis SAS
############################################################################
makeContrastType3SAS <- function(model, term, L)
{
  
  eps <- 1e-8
  #apply rule 1 (Goodnight 1976)
  
  #find all effects that contain term effect
  model.term <- terms(model)
  fac <- attr(model.term,"factors")
  names <- attr(model.term,"term.labels")
  classes.term <- attr(terms(model, FALSE), "dataClasses")#attr(model.term,"dataClasses")
  
  cols.eff <- which(colnames(L)==term)
  num.relate <- relatives(classes.term, term, names, fac)
  if( length(num.relate)==0 )
    colnums <- setdiff(1:ncol(L),cols.eff)
  if( length(num.relate)>0 )
  {
    cols.contain <- NULL
    for( i in 1:length(num.relate) )
      cols.contain <- c(cols.contain,which(colnames(L)==names[num.relate[i]]))
    colnums <- setdiff(1:ncol(L),c(cols.eff,cols.contain))   
  }
    
  for(colnum in colnums)
  {
    
    pivots <- which(abs(L[,colnum]) > eps)
    #pivots<-which(L[,colnum]!=0)
    if( length(pivots)>0 )
    {
      L[pivots[1],] <- L[pivots[1],]/L[pivots[1],colnum]
      nonzeros <- setdiff(pivots,pivots[1])
      if( length(nonzeros)!=0 )
      {
         for( nonzero in nonzeros )
         {
           L[nonzero,] <- L[nonzero,]-L[nonzero,colnum]*L[pivots[1],]
         }
      }
     
      L[pivots[1],] <- rep(0,ncol(L))
    }
  }
    
  nums <- which(apply(L,1,function(y) sum(abs(y)))!=0) 
  L <- L[nums,]
  
  if(is.vector(L))
    return(L)
  
  #orthogonalization
  if( length(cols.eff)>1 )
      zero.rows <- which(apply(L[,cols.eff],1,function(y) sum(abs(y)))==0)
  else
      zero.rows <- which(L[,cols.eff]==0)
      
  for(zero.row in zero.rows) 
  {
    w <- L[zero.row,]
    for(i in setdiff(1:nrow(L),zero.row))
    {
      if(sum(abs(L[i,]))!=0)
        L[i,] <- L[i,]-((w %*% L[i,])/(w %*% w)) %*% w
    }
    L[zero.row,] <- rep(0,ncol(L))
  }

  L[abs(L)<1e-6] <- 0
  
  nums <- which(apply(L,1,function(y) sum(abs(y)))!=0) 
  L <- L[nums,]
  return(L)
}

makeContrastType2 <- function(model, term, L, X.design, rho, fullCoefs){
  #find all effects that contain term effect
  model.term <- terms(model)
  fac <- attr(model.term,"factors")
  names <- attr(model.term,"term.labels")
  classes.term <- attr(terms(model, FALSE), "dataClasses")
  
  find.term <- which(colnames(X.design) == term)
  num.relate <- relatives(classes.term, term, names, fac)
  contain <- names[num.relate]
  if(length(contain) == 0 && (which(names == term) == length(names))){
    Lc <- L[find.term[which(find.term %in% rho$nums.Coefs)],]            
  }
  else{
    ind.indep <- which(colnames(X.design) != term & !(colnames(X.design) %in% contain))
    new.X <- cbind(X.design[,ind.indep], X.design[,-ind.indep])
    XtX <- crossprod(new.X)
    U <- doolittle(XtX)$U
    d <- diag(U)
    for(i in 1:nrow(U))
      if(d[i] > 0) U[i, ] <- U[i, ] / d[i]
    L <- U
    colnames(L) <- rownames(L) <-  
      c(names(fullCoefs)[ind.indep], names(fullCoefs)[-ind.indep])
    Lc <- L[which(colnames(new.X) == term), , drop = FALSE]
    Lc <- Lc[which(rownames(Lc) %in% names(rho$nums.Coefs)), , drop = FALSE]
    Lc <- Lc[, names(fullCoefs), drop = FALSE]
  }  
  return(Lc)
}

############################################################################
#get formula for model 
############################################################################
getFormula <- function(model, withRand=TRUE)
{
  fmodel <- formula(model)  
 
  if(withRand)
    return(fmodel)
  
  fm <- paste(fmodel)
  fmodel.red <- paste(fm[2],fm[1], 
                        paste(fm[3], 
                              paste(unlist(lapply(names(.fixedrand(model)$randeffs),
                                                  function(x) paste("(",x, ")"))), 
                                    collapse = " - "), 
                              sep = "-"))
  return(update(fmodel, fmodel.red))
}





#check if there are correlations between intercepts and slopes
checkCorr <- function(model)
{
   corr.intsl <- FALSE
   modelST <- getST(model)
   lnST <- length(modelST)
   for(i in 1:lnST)
   {    
      if(nrow(modelST[[i]])>1)
         corr.intsl <- TRUE
   } 
   return(corr.intsl) 
}




emptyAnovaLsmeansTAB <- function()
{
  result <- NULL
  anova.table <-  matrix(ncol=5,nrow=0)
  colnames(anova.table) <- c("Estimate","Standard Error", "DF", "F-value", "p-value")
  anova.table <- as.data.frame(anova.table)
  result$TAB.fixed <- anova.table
  lsmeans.summ <-  matrix(ncol=7,nrow=0)
  colnames(lsmeans.summ) <- c("Estimate","Standard Error", "DF", "t-value", "Lower CI", "Upper CI", "p-value")
  lsmeans.summ <- as.data.frame(lsmeans.summ)
  result$TAB.lsmeans <- lsmeans.summ
  return(result)
}



### check if there are no random terms in the model
checkPresRandTerms <- function(mf.final)
{
  sub.loc.rand <- substring.location(paste(mf.final)[3], "|")
  if(sub.loc.rand$first==0 && sub.loc.rand$last==0)
    return(FALSE)
  return(TRUE)
}

### compare mixed model versus fixed
compareMixVSFix <- function(model, mf.final, data, name.term)
{
 
  model.red <- refitLM(model)
  
  l.fix <- -2*logLik(model, REML=TRUE)[1]
  l.red <- -2*logLik(model.red, REML=TRUE)[1]
  
  p.chisq <- 1 - pchisq (l.red -l.fix ,1)
  infoForTerm <- saveInfoForTerm(name.term, l.red -l.fix, 1, p.chisq, 
                                 model.red = model.red)
  
  return(infoForTerm)
}



#create reduce slopes model
createModelRedSlopes <- function(x, term, fm, model, l.lmerTest.private.contrast)
{
  fm[3] <- paste(fm[3], "-", term, "+" , paste(x,collapse="+"))
  mf.final <-  as.formula(paste(fm[2],fm[1],fm[3], sep=""))
  mf.final <- update.formula(mf.final,mf.final)
  model.red <- updateModel(model, mf.final, getME(model, "is_REML"), 
                           l.lmerTest.private.contrast)
  return(model.red)
}



# ### get the random terms
getRandTerms <- function(fmodel)
{
  terms.fm <- attr(terms(fmodel),"term.labels")
  ind.rand.terms <- which(unlist(lapply(terms.fm,function(x) substring.location(x, "|")$first))!=0)
  return(unlist(lapply(terms.fm[ind.rand.terms],function(x) paste("(",x,")",sep=""))))
}




#####save results for fixed effects for model with only fixed effects
saveResultsFixModel <- function(result, model, type = 3)
{
  if(type==3)
    result$anova.table <- drop1(model, test="F")
  else{
    result$anova.table <- anova(model)
  }
  result$model <- model
  lsmeans.summ <-  matrix(ncol=7,nrow=0)
  colnames(lsmeans.summ) <- c("Estimate","Standard Error", "DF", "t-value", 
                              "Lower CI", "Upper CI", "p-value")
  result$lsmeans.table <- lsmeans.summ
  result$diffs.lsmeans.table <- lsmeans.summ
  return(result)
}



getREML <- function(model)
{
  if(inherits(model,"merMod"))
    return(getME(model, "is_REML"))

}

#update model
updateModel <- function(model, mf.final, reml.lmerTest.private, 
                        l.lmerTest.private.contrast, 
                        devFunOnly.lmerTest.private = FALSE)
{
  if(!mf.final == as.formula(paste(".~.")))
  {
    inds <-  names(l.lmerTest.private.contrast) %in% attr(terms(as.formula(mf.final)), 
                                                          "term.labels")
     #update contrast l.lmerTest.private.contrast
    l.lmerTest.private.contrast <- l.lmerTest.private.contrast[inds]
  }
 
 #data.update.lmerTest.private <- model.frame(model)  
 nfit <- update(object=model, formula.=mf.final, REML=reml.lmerTest.private ,
                contrasts=l.lmerTest.private.contrast, 
                devFunOnly = devFunOnly.lmerTest.private, 
 #               data = data.update.lmerTest.private,
                evaluate=FALSE)
 env <- environment(formula(model))
 assign("l.lmerTest.private.contrast", l.lmerTest.private.contrast, envir=env)
 assign("reml.lmerTest.private", reml.lmerTest.private, envir=env)
 assign("devFunOnly.lmerTest.private", devFunOnly.lmerTest.private, envir=env)
 #assign("data.update.lmerTest.private", data.update.lmerTest.private, envir=env)
 nfit <- eval(nfit, envir = env) 
 return(nfit)   
}


### Rhune's code for making type 1 SS
doolittle <- function(x, eps = 1e-6) {
  
  if(!is.matrix(x)) stop("argument 'x' is not a matrix")
  if(ncol(x) != nrow(x))
    stop( "argument x is not a square matrix" )
  if (!is.numeric(x) )
    stop( "argument x is not numeric" )
  n <- nrow(x)
  L <- U <- matrix(0, nrow=n, ncol=n)
  diag(L) <- rep(1, n)
  for(i in 1:n) {
    ip1 <- i + 1
    im1 <- i - 1
    for(j in 1:n) {
      U[i,j] <- x[i,j]
      if (im1 > 0) {
        for(k in 1:im1) {
          U[i,j] <- U[i,j] - L[i,k] * U[k,j]
        }
      }
    }
    if ( ip1 <= n ) {
      for ( j in ip1:n ) {
        L[j,i] <- x[j,i]
        if ( im1 > 0 ) {
          for ( k in 1:im1 ) {
            L[j,i] <- L[j,i] - L[j,k] * U[k,i]
          }
        }
        L[j, i] <- if(abs(U[i, i]) < eps) 0 else L[j,i] / U[i,i]
      }
    }
  }
  L[abs(L) < eps] <- 0
  U[abs(U) < eps] <- 0
  list( L=L, U=U )
}


refitLM <- function(obj, l.lmerTest.private.contrast="contr.SAS") {

  mm <- model.frame(obj)
  colnames(mm)[1] <- "y"
  fo <- getFormula(obj, withRand=FALSE)# formula(obj,fixed.only=TRUE)
  if(fo != as.formula(.~.))
  {
    inds <-  names(l.lmerTest.private.contrast) %in% attr(terms(fo), "term.labels")
    #update contrast l
    l.lmerTest.private.contrast <- l.lmerTest.private.contrast[inds]
  }
  fo <- update(fo, y ~ .)
  lm(fo, data=mm, contrasts = l.lmerTest.private.contrast)
}

###########################################################
# fill anova table
###########################################################
fillAnovaTable <- function(result, anova.table)
{
  for (i in 1:length(result))
  {
    if(!result[[i]]$name %in% rownames(anova.table))
      next
    anova.table[result[[i]]$name, 4] <- result[[i]]$denom
    anova.table[result[[i]]$name, 5] <- result[[i]]$Fstat
    anova.table[result[[i]]$name, which(colnames(anova.table)=="Pr(>F)")] <- 
      result[[i]]$pvalue
    if(!is.na(result[[i]]$ss)){
      anova.table[result[[i]]$name, "Sum Sq"] <- result[[i]]$ss
      anova.table[result[[i]]$name, "Mean Sq"] <- result[[i]]$ms
    }
    
    
    
  }
  anova.table
}



# format table according to elim.num column
formatElimNumTable <- function(table)
{
  if("elim.num" %in% colnames(table)){
    table[which(table[,"elim.num"]==0),"elim.num"] <- 1000
    table <- table[with(table, order(elim.num, decreasing=FALSE)),]
    table[,"elim.num"] <- as.character(table[,"elim.num"])
    table[which(table[,"elim.num"]=="1000"),"elim.num"] <- "kept"
  }
  return(table)
}


updateAnovaTable <- function(resNSelim){
  anm <- anova(resNSelim$model)
  anova.table <- resNSelim$anova.table
  anova.table[rownames(anm), c("Sum Sq", "Mean Sq", "NumDF")] <-
    as.matrix(anm[, c("Sum Sq", "Mean Sq", "Df")])
  anova.table
}


.getKeepInter <- function(model.eff, split.keep.eff){
  if(setequal(unlist(strsplit(model.eff, ":")), split.keep.eff))
    return(model.eff)
}

.findKeepEff <- function(keep.eff, model.effs){
 
  ## for main terms
  if(length(grep(":", keep.eff)) == 0){
    if(keep.eff %in% model.effs)
      return(keep.eff)
  }
  else{
    ## for interaction terms
    return(unlist(lapply(model.effs, .getKeepInter, 
                         unlist(strsplit(keep.eff, ":")))))
  }   
}

.findKeepEffSlope <- function(keep.eff, randTermsSl){
  
  ind.Slope <- unlist(lapply(randTermsSl, function(x) x$sl.part == keep.eff)) 
  if(sum(ind.Slope) == 0){
    ind.Gr <- unlist(lapply(randTermsSl, function(x) x$gr.part == keep.eff)) 
    if(sum(ind.Gr) > 0)
      return(c(names(randTermsSl[ind.Slope]), keep.eff))
  } 
  c(names(randTermsSl[ind.Slope]))
}

.findKeepRandEff <- function(keep.eff, model.effs){
  randTerms <- getRandTermsTable(names(model.effs))
  is.scalar <- unlist(lapply(randTerms, function(x) x$sl.part == "1"))
  randTermsScal <- randTerms[is.scalar]
  randTermsSl <- randTerms[!is.scalar]
  keep.eff.scal <- .findKeepEff(keep.eff, names(randTermsScal))
  keep.eff.slope <- .findKeepEffSlope(keep.eff, randTermsSl)
  c(keep.eff.scal, keep.eff.slope)  
}

.getKeepEffs <- function(keep.effs, model.effs){
  #rand.terms.table <- getRandTermsTable(rand.terms)
  randeffs <- unlist(lapply(keep.effs, .findKeepRandEff, 
                            unlist(model.effs$randeffs)))
  fixedeffs <- unlist(lapply(keep.effs, .findKeepEff, 
                             unlist(model.effs$fixedeffs)))
  return(list(randeffs = randeffs, fixedeffs = fixedeffs))
}


## the same function is in SensMixed package
## list of random and fixed terms of the model
.fixedrand <- function(model)
{  
  effs <- attr(terms(formula(model)), "term.labels")
  neffs <- length(effs)
  randeffs <- effs[grep(" | ", effs)]
  randeffs <- sapply(randeffs, function(x) substring(x, 5, nchar(x)))
  fixedeffs <- effs[!(effs %in% names(randeffs))]
  return(list(randeffs=randeffs, fixedeffs=fixedeffs))
}

