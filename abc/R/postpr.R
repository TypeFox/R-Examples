######################################################################
#
# postpr.R
#
# copyright (c) 2011-05-30, Katalin Csillery, Olivier Francois and
# Michael GB Blum with some initial code from Mark Beaumont
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
# 
# Part of the R/abc package
# Contains: postpr, summary.postpr
#
######################################################################

postpr <- function(target, index, sumstat, tol, subset=NULL, method, corr=TRUE, kernel="epanechnikov",
                   numnet = 10, sizenet = 5, lambda = c(0.0001,0.001,0.01), trace = TRUE, maxit = 500, ...){

  linout <- FALSE
  call <- match.call()

  ## general checks that the function is used correctly
  ## ###################################################

  if(missing(target)) stop("'target' is missing with no default", call.=F)
  if(missing(index)) stop("'index' is missing with no default", call.=F)
  if(missing(sumstat)) stop("'sumstat' is missing with no default", call.=F)
  if(!is.vector(index)) stop("'index' has to be a vector.", call.=F)
  if(!is.matrix(sumstat) && !is.data.frame(sumstat) && !is.vector(sumstat)) stop("'sumstat' has to be a matrix, data.frame or vector.", call.=F)
  if(missing(tol)) stop("'tol' is missing with no default", call.=F)
  if(missing(method)) stop("'method' is missing with no default", call.=F)
  if(!any(method == c("rejection", "mnlogistic", "neuralnet")))
    stop("Method must be 'rejection', 'mnlogistic' or 'neuralnet'.", call.=F)
  if(length(unique(index)) == 1)
    stop("At least two different models must be given.", call.=F)
  if(method == "rejection") rejmethod <- TRUE
  else rejmethod <- FALSE

  if(is.data.frame(sumstat)) sumstat <- as.matrix(sumstat)
  if(is.vector(sumstat)) sumstat <- matrix(sumstat, ncol=1)
  if(is.list(target)) target <- unlist(target)
  if(length(target)!=dim(sumstat)[2]) stop("Number of summary statistics in 'target' has to be the same as in 'sumstat'.", call.=F)
  if(length(index) != length(sumstat[,1]))
    stop("'index' must be the same length as the number of rows in 'sumstat'.", call.=F)
  
  if(!any(kernel == c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine"))){
    kernel <- "epanechnikov"
    warning("Kernel is incorrectly defined. Setting to default kernel (Epanechnikov)", call.=F, immediate=T)
  }
  
  if(is.vector(sumstat)) sumstat <- matrix(sumstat, ncol=1)
  if(length(target)!=dim(sumstat)[2]) stop("Number of summary statistics in 'target' has to be the same as in 'sumstat'.", call.=F)

  index <- factor(index)
  mymodels <- levels(index)
  if(!is.numeric(target)) target <- as.numeric(target)
  
  ## parameter and/or sumstat values that are to be excluded
  ## #######################################################
  gwt <- rep(TRUE,length(sumstat[,1]))
  gwt[attributes(na.omit(sumstat))$na.action] <- FALSE
  if(is.null(subset)) subset <- rep(TRUE,length(sumstat[,1]))
  gwt <- as.logical(gwt*subset)
  
  sumstat <- as.data.frame(sumstat) ## ?????????????????
  ## extract names of statistics if given
  ## ####################################
  nss <- length(sumstat[1,])
  if(!length(colnames(sumstat))){
    warning("No summary statistics names are given, using S1, S2, ...", call.=F, immediate=T)
    statnames <- paste("S", 1:nss, sep="")
  }
  else statnames <- colnames(sumstat)
    
  ## stop if zero var in sumstat
  ## ###########################
  cond1 <- as.logical(apply(sumstat, 2, function(x) length(unique(x))-1))
  if(!all(cond1)) stop("Summary statistic(s) have zero variance.", call.=F)
  if(!any(cond1)){
    warning("Statistic(s) ", statnames[!cond1], " have zero variance. Excluding from estimation....", sep="\t", call.=F, immediate=T)
    sumstat <- sumstat[,cond1]
    nss <- length(sumstat[1,])
    statnames <- colnames(sumstat)
    target <- target[cond1]
  }

  ## scale everything
  ## ################
  scaled.sumstat <- sumstat
  for(j in 1:nss){
    scaled.sumstat[,j] <- normalise(sumstat[,j],sumstat[,j][gwt])
  }
  
  for(j in 1:nss){
    target[j] <- normalise(target[j],sumstat[,j][gwt])
  }

  ## calculate euclidean distance
  ## ############################
  sum1 <- 0
  for(j in 1:nss){
    sum1 <- sum1 + (scaled.sumstat[,j]-target[j])^2
  }
  dist <- sqrt(sum1)

  # includes the effect of gwt in the tolerance
  dist[!gwt] <- floor(max(dist[gwt])+10)

  # wt1 defines the region we're interested in
  abstol <- quantile(dist,tol)
  if(kernel == "gaussian") wt1 <- rep(TRUE, length(dist)) ## ???????????????????????????
  else{
      ceiling(length(dist)*tol)->nacc
      sort(dist)[nacc]->ds
      wt1 <- (dist <= ds)
      aux<-cumsum(wt1)
      wt1 <- wt1 & (aux<=nacc)
  }
  
  ## select summary statistics in region
  ## ##################################
  ss <- scaled.sumstat[wt1,]
  values <- index[wt1]
  pred <- table(values)/length(values)

  statvar <- as.logical(apply(scaled.sumstat[wt1, , drop=FALSE], 2, function(x) length(unique(x))-1))
  cond2 <- !any(statvar)
  
  if(cond2 && !rejmethod)
    stop("Zero variance in the summary statistics in the selected region.\nTry: checking summary statistics, choosing larger tolerance, or rejection method.", call.=F)

  
  ## if simple rejection or in the selected region there is no var in sumstat
  ## #########################################################################
  if(rejmethod){
    if(cond2) warning("Zero variance in the summary statistics in the selected region. Check summary statistics, consider larger tolerance.", call.=F, immediate=T)
    weights <- NULL
    pred.logit <- NULL
  }

  ## regression correction
  ## ######################
  else{
    if(cond2) cat("Warning messages:\nStatistic(s)",
                  statnames[!statvar],
                  "has/have zero variance in the selected region.\nConsider using larger tolerance or the rejection method or discard this/these statistics.\n", sep="\t")
    
    ## weights
    if(kernel == "epanechnikov") weights <- 1 - (dist[wt1]/abstol)^2
    if(kernel == "rectangular") weights <- dist[wt1]/abstol
    if(kernel == "gaussian") weights <- 1/sqrt(2*pi)*exp(-0.5*(dist/abstol)^2)
    if(kernel == "triangular") weights <- 1 - abs(dist[wt1]/abstol)
    if(kernel == "biweight") weights <- (1 - (dist[wt1]/abstol)^2)^2
    if(kernel == "cosine") weights <- cos(pi/2*dist[wt1]/abstol)

    ok <- index[wt1] # models accepted
    fml <- as.formula(paste("ok ~ ", paste(statnames, collapse= "+")))
      
    if(length(unique(ok)) < length(mymodels)) {
      warning(paste("There are",length(mymodels),"models but only",length(unique(ok)), "for which simulations have been accepted.\nNo regression is performed, method is set to rejection.\nConsider increasing the tolerance rate."),sep="", call.=F, immediate=T)
      weights <- NULL
      pred.logit <- NULL
      method <- "rejection"
    }

    ## calculating the number of weights for multinom
    mymnw <- (nss+2) * length(mymodels)
    
    if(method == "mnlogistic"){
        ss<-data.frame(ss)
        colnames(ss)<-statnames
      fit1 <- multinom(fml, data = ss, weigths = weights, trace=F, MaxNWts = mymnw + 1, ...)
      target <- as.data.frame(matrix(target, nrow=1))
      names(target) <- statnames
      pred <- predict(fit1, target, type="probs")
      if(length(pred) == 1){
        pred <- c(1-pred,pred)
        names(pred) <- levels(ok)
      }
    }
    
    else if(method == "neuralnet"){
        ss<-data.frame(ss)
        colnames(ss)<-statnames
      lambda <- sample(lambda, numnet, replace=T)
      target <- as.data.frame(matrix(target, nrow=1))
      names(target) <- statnames
      pred <- 0
      for(i in 1:numnet){
        fit1 <- nnet(fml, data=ss, weights = weights, decay = lambda[i],
                     size = sizenet, trace = trace, linout = linout, maxit = maxit, ...)
        if(length(mymodels)==2)
          {
            auxm<-predict(fit1, target, type="raw")
            pred <- pred + c(1-auxm,auxm)
          }
        else
          pred <- pred + predict(fit1, target, type="raw")
      }
      pred <- pred/numnet
      if(length(mymodels)!=2)
        {
          temp <- rep(0, length(mymodels))
          names(temp) <- mymodels
          temp[match(colnames(pred), mymodels)] <- pred
          pred <- temp
        }
      else names(pred) <- levels(ok)
    }
    ## correction for potentially different numbers of simulations per models
    ratio <- (pred*length(index)*tol) / table(index)
    pred <- ratio/sum(ratio)
    attributes(dimnames(pred)) <- NULL
  }
  
  if(rejmethod){
    postpr.out <- list(values=values, ss=ss, call=call, na.action=gwt, method=method, corr=corr, nmodels=table(index),
                       numstat=nss, names=list(models=mymodels, statistics.names=statnames))
  }
  else{
    postpr.out <- list(values=values, pred=pred, ss=ss, weights=weights, call=call, na.action=gwt, method=method, corr=corr, nmodels=c(table(index)),
                       numstat=nss, names=list(models=mymodels, statistics.names=statnames))
  }
  class(postpr.out) <- "postpr"
  invisible(postpr.out)
}

summary.postpr <- function(object, rejection = TRUE, print = TRUE, digits = max(3, getOption("digits")-3), ...){

  if (!inherits(object, "postpr")) 
    stop("Use only with objects of class \"postpr\".", call.=F)
  
  postpr.out <- object
  cl <- postpr.out$call
  npost <- length(postpr.out$values)
  pred <- postpr.out$pred
  allvals <- postpr.out$values
  postmod <- levels(postpr.out$values)
  nmod <- length(postmod)
  method <- postpr.out$method
  corr <- postpr.out$corr
  nmodels <- postpr.out$nmodels
  
  if(print){
      cat("Call: \n")
      dput(cl, control=NULL)  
      cat(paste("Data:\n postpr.out$values (",npost," posterior samples)\n", sep=""))
      cat(paste("Models a priori:\n "))
      cat(postpr.out$names$models, sep=", ")
      cat(paste("\nModels a posteriori:\n "))
      cat(postmod, sep=", ")
      if(corr & length(unique(nmodels))>1){
          cat("\n")
          warning("Posterior model probabilities are corrected for unequal number of simulations per models.", immediate.=T, call.=F)
      }
      cat("\n\n")
  }
  
  if(rejection || method == "rejection"){
      
      if(print) cat("Proportion of accepted simulations (rejection):\n")
      allpr <- table(allvals)/length(allvals)
      if(corr){
          ratio <- (allpr*npost) / nmodels
          allpr <- ratio/sum(ratio)
      }
      prnames <- dimnames(allpr)$allvals
      allpr <- c(allpr); names(allpr) <- prnames
      if(print) print(round(allpr, digits=digits))
      
      if(nmod>1){
          pr.rej <- table(allvals)/length(allvals)
          bf.rej <- t(matrix(pr.rej, nmod, nmod, byrow=T)/matrix(pr.rej, nmod, nmod, byrow=F))
          colnames(bf.rej) <- postmod
          rownames(bf.rej) <- postmod
          bf.rej <- as.table(bf.rej)
          if(print){
              cat("\nBayes factors:\n")
              print(round(bf.rej, digits=digits))
              cat("\n\n")
          }
      }
      else bf.rej <- NA
      
  }
  
  if(method == "mnlogistic" | method == "neuralnet"){
    
    if(print){
      cat(paste("Posterior model probabilities (", method, "):\n", sep=""))
      print(round(pred, digits=digits))
    }
    if(nmod>1){
      bf.reg <- t(matrix(pred[pred!=0], nmod, nmod, byrow=T)/matrix(pred[pred!=0], nmod, nmod, byrow=F))
      colnames(bf.reg) <- postmod
      rownames(bf.reg) <- postmod
      bf.reg <- as.table(bf.reg)
      if(print){
        cat("\nBayes factors:\n")
        print(round(bf.reg, digits=digits))
        cat("\n")
      }
    }
    else bf.reg <- NA
      
    if(rejection){
      if(method == "mnlogistic")
        out <- list(rejection=list(Prob=allpr, BayesF=bf.rej), mnlogistic=list(Prob=pred, BayesF=bf.reg))
      if(method == "neuralnet")
        out <- list(rejection=list(Prob=allpr, BayesF=bf.rej), neuralnet=list(Prob=pred, BayesF=bf.reg))
    }
    else{
      out <- list(Prob=pred, BayesF=bf.reg)
    }
  }
  else out <- list(Prob=allpr, BayesF=bf.rej)
  invisible(out)
}
  


