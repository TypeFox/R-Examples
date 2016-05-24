#  File R/ergm.ego.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2016 Statnet Commons
#######################################################################
ergm.ego <- function(formula, popsize=1, offset.coef=NULL, ..., control=control.ergm.ego(), na.action=na.fail, do.fit=TRUE){
  check.control.class()
  
  stats.est <- control$stats.est
  stats.wt <- control$stats.wt
  egodata <- get(as.character(formula[[2]]), envir=environment(formula))

  sampsize <- dim(egodata)[1]
  ppopsize <-
    if(is.numeric(control$ppopsize)) control$ppopsize
    else switch(control$ppopsize,
                auto = if(missing(popsize) || popsize==1) sampsize*control$ppopsize.mul else popsize*control$ppopsize.mul,  
                samp = sampsize*control$ppopsize.mul,
                pop = popsize*control$ppopsize.mul)
  
  if(ppopsize < sampsize) stop("Using a smaller pseudopopulation size than sample size does not make sense.")
  else if(ppopsize == sampsize && !is.null(egodata$egoWt) && var(egodata$egoWt)>sqrt(.Machine$double.eps))
    warning("Using pseudopoulation size equal to sample size under weighted sampling: results may be highly biased. Recommend increasing popsize.mul control parameter.")
  
  message("Constructing pseudopopulation network.")
  popnw <- as.network(egodata, ppopsize, scaling=control$ppop.wt)
  if(network.size(popnw)!=ppopsize){
    message("Note: Constructed network has size ", network.size(popnw), ", different from requested ", ppopsize,". Estimation should not be meaningfully affected.")
    ppopsize <- network.size(popnw)
  }

  w <- switch(stats.wt,
              data=egodata$egoWt,
              ppop=tabulate(popnw %v% "ego.ind", nbins=nrow(egodata))
              )
  
  # Get the sample h values.
  stats <- try(summary(remove.offset.formula(formula), individual=TRUE))
  if(!inherits(stats,"try-error")){
    # h is just a matrix, so this will do the sensible thing.
    tmp <- na.action(cbind(w,stats))
    w <- tmp[,1]
    stats <- tmp[,-1, drop=FALSE] * ppopsize
    n <- length(w)
    
    wmean <- function(w,s){
      w <- w/sum(w)
      colSums(s*w)
    }
    m <- wmean(w,stats)
    
    if(stats.est=="bootstrap"){
      m.b <- t(replicate(control$boot.R,{
                           i <- sample.int(length(w),replace=TRUE)
                           wmean(w[i],stats[i,,drop=FALSE])
                         }))
      m <- m - (colMeans(m.b)-m)
      
    }else if(stats.est=="jackknife"){
      m.j <- t(sapply(seq_len(n), function(i){
                        wmean(w[-i],stats[-i,,drop=FALSE])
                      }))
      m <- n*m - (n-1)*colMeans(m.j)
    }
    
    # TODO: Include finite-population correction here:
    v <- switch(stats.est,
                bootstrap = cov(m.b),
                jackknife = (n-1)/n*crossprod(sweep(m.j,2,colMeans(m.j))),
                naive = {w <- w/sum(w); crossprod(sweep(stats, 2, m, "-")*sqrt(w))/(1-sum(w^2))*sum(w^2)},
                asymptotic = .asymptotic.var(stats, w)/length(w)
                )
  }else{
    if(stats.est %in% c("naive","asymptotic"))
      stop("Non-scaling statistic detected: use bootstrap or jackknife variance estimator.")
    if(do.fit && popsize!=ppopsize)
      warning("Non-scaling statistic detected when trying to fit a model: network-size invariant parametrization probably does not exist so pseudopopulation size should equal the population size.")

    n <- nrow(egodata)
    m <- summary(remove.offset.formula(formula), basis=egodata, individual=FALSE, scaleto=ppopsize)
      
    if(stats.est=="bootstrap"){
      m.b <- t(replicate(control$boot.R,{
                           i <- sample.int(length(w),replace=TRUE)
                           e <- egodata[i,]
                           summary(remove.offset.formula(formula), basis=e, individual=FALSE, scaleto=ppopsize)
                         }))
      m <- m - (colMeans(m.b)-m)
      
    }else if(stats.est=="jackknife"){
      m.j <- t(sapply(seq_len(n), function(i){
                        e <- egodata[-i,]
                        summary(remove.offset.formula(formula), basis=e, individual=FALSE, scaleto=ppopsize)
                      }))
      m <- n*m - (n-1)*colMeans(m.j)
    }
    
    # TODO: Include finite-population correction here:
    v <- switch(stats.est,
                bootstrap = cov(m.b),
                jackknife = (n-1)/n*crossprod(sweep(m.j,2,colMeans(m.j)))
                )
  }
  
  ergm.formula <- ergm.update.formula(formula,popnw~offset(netsize.adj)+.,from.new="popnw")
  ergm.offset.coef <- c(-log(ppopsize/popsize),offset.coef)
  out <- list(v=v, m=m, formula=formula, ergm.formula=ergm.formula, offset.coef=offset.coef, ergm.offset.coef=ergm.offset.coef, egodata=egodata, ppopsize=ppopsize, popsize=popsize)
  
  if(do.fit){

    ergm.fit <- ergm(ergm.formula, target.stats=m, offset.coef=ergm.offset.coef,..., eval.loglik=FALSE,control=control$ergm.control)

    ## Workaround to keep mcmc.diagnostics from failing. Should be removed after fix is released.
    if(inherits(ergm.fit$sample,"mcmc.list")){
      for(thread in 1:nchain(ergm.fit$sample))
        ergm.fit$sample[[thread]][,1] <- 0
    } else ergm.fit$sample[,1] <- 0
    ergm.fit$drop[1] <- 0

    coef <- coef(ergm.fit)

    oi <- offset.info.formula(ergm.formula)

    DtDe <- -ergm.fit$hessian[!oi$theta,!oi$theta,drop=FALSE]

    vcov <- matrix(NA, length(coef), length(coef))
  
    vcov[!oi$theta,!oi$theta] <- solve(DtDe)%*%v%*%solve(DtDe)
    
    rownames(vcov) <- colnames(vcov) <- names(coef)

    out <- c(out, list(covar=vcov, ergm.covar=ergm.fit$covar, DtDe=DtDe))
    out <- modifyList(ergm.fit, out)
  }
  class(out) <- c("ergm.ego","ergm")
  out
}

.asymptotic.var <- function(X, w){
  n <- length(w)
  p <- ncol(X)
  wX <- cbind(w,sweep(X,1,w,"*"))
  m <- colSums(wX[,-1,drop=FALSE])/sum(w)
  S <- cov(wX)
  A <- 1/mean(w)*cbind(-m,diag(1,nrow=p))
  A%*%S%*%t(A)
}

