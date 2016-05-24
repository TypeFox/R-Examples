#  file likelihoodAsy.R  (various functions)
#  This file is a component of the package 'likelihoodAsy' for R 
#  copyright (C) 2014 Ruggero Bellio and Donald A. Pierce.
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#-----------------------------------------------------------------

.newscores <- function(p, k, C.hat, C.til, score.hat, score.til)
 { 
#
# Internal function which transforms the theta-scores to those for (psi, lam)
# with lam being a suitable p-1 coordinates of theta 
# 
    A.hat <- C.hat[k]   ## part psi / part theta_k
    B.hat <- C.hat[-k]  ## part psi / part lambda
    A.til <- C.til[k]   ## part psi / part theta_k
    B.til <- C.til[-k]  ## part psi / part lambda
    score.new.hat <- score.new.til <- rep(0,p)
    score.new.hat[k] <- score.hat[k] / A.hat
    score.new.hat[-k] <- score.hat[k]*(-1)*B.hat/A.hat + score.hat[-k]
    score.new.til[k] <- score.til[k] / A.til
    score.new.til[-k] <- score.til[k]*(-1)*B.til/A.til + score.til[-k]
    ans <- list(score.new.hat=score.new.hat, score.new.til=score.new.til)
    return(ans)
}



.newinfo <- function(p, k, C.hat, C.til, j.hat, j.til, score.hat.data, score.til.data, theta.hat, theta.til, fpsi)
{
#
#  Internal function which transforms from obs info for theta to that for (psi,lam) 
#  with lam being a suitable p-1 coordinates of theta 
#
    pm1 <- p-1
	A.hat <- C.hat[k]   ## part psi / part theta_k
    B.hat <- C.hat[-k]  ## part psi / part lambda
    A.til <- C.til[k]   ## part psi / part theta_k
    B.til <- C.til[-k]  ## part psi / part lambda
    ## for the "hat" fit this is standard transformation stuff
	dnew.dold <- matrix(0, nrow=p,ncol=p)
	dnew.dold[k,k] <- matrix(A.hat,ncol=1,nrow=1)
	dnew.dold[k,-k] <- matrix(rep(0,pm1),nrow=1,ncol=pm1)
	dnew.dold[-k,k] <- matrix(B.hat,nrow=pm1,ncol=1)
	dnew.dold[-k,-k] <- diag(pm1)
	T <- solve(dnew.dold) 
	j.hat.new <- T %*% j.hat %*% t(T) 
    ## special calc for observed info at tilde
	## this added term arises because the score is not zero at tilde
    ##  need to transform to score for (psi, lam)
	score.new.til.data.k <- score.til.data[k] / A.til   ##right side uses old param theta
	obj.score <- .newscores(p,k,C.hat,C.til, score.hat=score.hat.data,score.til.data)
	score.new.til.data <- obj.score$score.new.til 
    ## get the new T at tilde
	dnew.dold <- matrix(0, nrow=p,ncol=p)
	dnew.dold[k,k] <- matrix(A.til,ncol=1,nrow=1)
	dnew.dold[k,-k] <- matrix(rep(0,pm1),nrow=1,ncol=pm1)
	dnew.dold[-k,k] <- matrix(B.til,nrow=pm1,ncol=1)
	dnew.dold[-k,-k] <- diag(pm1)
	T_til <- solve(dnew.dold) 
    added.term <- score.new.til.data[k] * pracma::hessian(fpsi,theta.til)   
	j.til.new <- T_til %*% (j.til + added.term) %*% t(T_til) 
	ans <- list(j.hat.new = j.hat.new, j.til.new = j.til.new)
    return(ans)
	}



rstar <- function(data, thetainit, floglik, fscore=NULL, fpsi, psival, datagen, R=1000, seed=NULL, trace=TRUE, ronly=FALSE, psidesc=NULL, constr.opt="solnp")
{
#
#  r* computation
#  
  ####### checking input arguments
  if(!is.list(data))
    {
     warning("data should be provided as a list\n")
     data <- as.list(data)
    } 
  if(!is.numeric(thetainit))
    stop("a starting point for the parameter theta is required \n")
  if(!is.numeric(psival))
    stop("a value for the parameter of interest is required \n")
  f0 <- floglik(thetainit, data)
  if(!is.numeric(f0))
     stop("problems in the loglikelihood function \n")
  if(!is.null(fscore))
   {
   	g0 <- fscore(thetainit, data)
    if(length(g0)!=length(thetainit))
      stop("size of starting point different from the size of score function \n")
   }
  if(!ronly)
    {
     data0 <- datagen(thetainit, data)
     f0 <- floglik(thetainit, data0)
     if(!is.numeric(f0))
     stop("problems in the function to simulate data \n")
    } 
   if( (constr.opt!="solnp") & (constr.opt!="alabama"))
   stop("constrained optimizer must be either 'solnp' or 'alabama'") 
    
  ####### getting the MLE
  p <- length(thetainit)
  if(trace) cat("get mle ....","\t")
  min.floglik <- function(theta, data) (-1) * floglik(theta,data)
  min.fscore <- if(is.null(fscore)) NULL else function(theta, data) (-1) * fscore(theta,data)
  obj.hat <- nlminb(thetainit, min.floglik, min.fscore, data=data )
  theta.hat <- obj.hat$par
  el.hat <- floglik(theta.hat, data)
  if(!ronly) 
     { j.hat <- if(is.null(fscore)) -pracma::hessian(floglik, theta.hat, data=data)
                else  -pracma::jacobian(fscore, theta.hat, data=data)
       score.hat.data <- if(is.null(fscore)) pracma::grad(floglik, theta.hat, data=data)
	                     else fscore(theta.hat, data=data)
     }
  ##### p>1: with NPs 
  if(length(thetainit)>1)
   {
   	if(!ronly)
	   {
	   	var.theta.hat <- solve(j.hat)
        se.theta.hat <- sqrt(diag(var.theta.hat))
	   }
   	if(trace) cat("get mle under the null....","\n")
    psifcn.mod <- if (constr.opt=="solnp") function(theta, data) fpsi(theta)  	
                  else function(theta, data) fpsi(theta) - psival 
  	objHyp <- if (constr.opt=="solnp") solnp(theta.hat,fun=min.floglik, eqfun=psifcn.mod, eqB=psival,
	           control=list(trace=0), data=data)
	          else constrOptim.nl(theta.hat,fn=min.floglik, heq=psifcn.mod, gr=min.fscore,
	           control.outer=list(trace=FALSE), data=data) 
    theta.til <- objHyp$par  
    el.til <- floglik(theta.til,data)
   	psi.hat <- fpsi(theta.hat)
   	if(!ronly)
	   {
        j.til <-  if(is.null(fscore)) -pracma::hessian(floglik, theta.til, data=data) 
                  else -pracma::jacobian(fscore, theta.til,data=data)
        score.til.data <- if(is.null(fscore)) pracma::grad(floglik,theta.til,data=data)
	                      else fscore(theta.til,data)
	    ## other main results from fit
        dpsi.dtheta <- pracma::grad(fpsi, theta.hat) 
        var.psi.hat <- dpsi.dtheta %*% var.theta.hat %*% dpsi.dtheta
        se.psi.hat <- sqrt(var.psi.hat)             
	   }
	 r <-  sqrt(2 * (el.hat - el.til)) * sign(psi.hat - psival)
	 ## find p-1 coords to use as NP
	 ## for NP omit any coord theta with nonzero grad element
	 if(!ronly)
	   {
	    C.hat <- pracma::grad(fpsi, theta.hat)
	    C.til <- pracma::grad(fpsi, theta.til)
        k <- which(C.hat != 0)[1]   # the index of first nonzero element  
       obj.info <- .newinfo(p, k, C.hat, C.til, j.hat, j.til, score.hat.data, score.til.data, theta.hat, theta.til, fpsi)
	   j.hat.new <- obj.info$j.hat.new
	   j.til.new <- obj.info$j.til.new 
	   }
    if(ronly) out<- list(r=r, theta.hat=theta.hat,  psi.hat = psi.hat, theta.hyp=theta.til, psi.hyp=psival)
    if(!ronly) 
     {
      ####### perform simulations 
    	if(trace) cat("start Monte Carlo computation","\n")
	    meanAll <- rep(0, 2*p+1) 
	    prodAll <- matrix(0, 2*p+1, 2*p+1)
	    dataSim <- as.list(data,all.names=TRUE)
	    seed.in <- if(is.null(seed)) sample.int(2^30,1)
                   else seed  
	    set.seed(seed.in)
	    if(trace) pb <- txtProgressBar(style=3) 
	    for (i in 1:R)
         { 
           if( (i%%(R/10)==0) & trace) setTxtProgressBar(pb, i/R)
           dataSim <- datagen(theta.hat, data=data)  
           l1 <- floglik(theta.hat, dataSim)
           l0 <- floglik(theta.til, dataSim)
	       score.hat <-  if(is.null(fscore))  pracma::grad(f=floglik, x0=theta.hat, data=dataSim)    
	                     else fscore(theta.hat, dataSim)
	       score.til <-  if(is.null(fscore)) pracma::grad(f=floglik, x0=theta.til, data=dataSim) 
                          else fscore(theta.til, dataSim)
	       ## get (psi,lam) scores from theta scores
           obj.score <- .newscores(p, k, C.hat, C.til, score.hat, score.til)
	       uhh <- c(obj.score$score.new.hat,  obj.score$score.new.til, (l1-l0))    
           meanAll <- meanAll + uhh / R
           prodAll <- prodAll +  tcrossprod(uhh) / R             
          }
      if(trace) close(pb)     
      #######  computes SS approximations 
      covAll <- prodAll * R/(R-1) - tcrossprod(meanAll)  * R / (R-1)
      S <-  covAll[1:p, (p+1):(2*p)] 
      i.hat <- covAll[1:p, 1:p]       
      i.hatInv <- qr.solve(i.hat,tol=10^-20)
      q <-  covAll[2*p+1, 1:p]        
      SS2 <- t(S) %*% i.hatInv %*% j.hat.new 
      SS1 <- q %*% i.hatInv %*% j.hat.new
      #### standard r* computation 
      indpsi <- k  ##  this is now the coord of theta which is omitted to leave lambda
      numU <- SS1[indpsi] -  SS2[-indpsi,indpsi] %*% solve(t(SS2[-indpsi,-indpsi])) %*% SS1[-indpsi] 
      j.hatInv.new <- solve(j.hat.new)  
      jProf.new <- 1 / j.hatInv.new[indpsi, indpsi] 
      u <- numU / sqrt(jProf.new)
      CPsi <- det(as.matrix(SS2[-indpsi, -indpsi])) / sqrt(det(as.matrix(j.til.new[-indpsi,-indpsi]   )) * det(as.matrix(j.hat.new[-indpsi, -indpsi]))) 
      NP <- (1/r) * log(CPsi)  
      INF <- (1/r) * log(u/r)   
      rs <- r + NP + INF 
      out <- list(r=r, NP=drop(NP), INF=drop(INF), rs=drop(rs), theta.hat = theta.hat, info.hat=j.hat, se.theta.hat=se.theta.hat,
                   psi.hat = psi.hat, se.psi.hat = drop(se.psi.hat), theta.hyp=theta.til, psi.hyp=psival, seed=seed.in)
     }
  } 
  ##### p=1: no NPs 
  if(length(thetainit)==1)
   {
   	if(!ronly)
   	  {
   	   j.hat <- as.numeric(j.hat)
   	   var.theta.hat <- 1/(j.hat)
       se.theta.hat <- sqrt(var.theta.hat)
   	  } 
   	psifcn.mod <- function(theta, data) fpsi(theta) - psival   
   	#### solve for theta=fpsi(psival)^-1	
  	objHyp <- nleqslv(theta.hat, psifcn.mod, data=data)
	theta.til <- objHyp$x
    el.til <- floglik(theta.til,data) 
    psi.hat <- fpsi(theta.hat)  
    if(!ronly)
   	  { 
       j.til <-  if(is.null(fscore)) -pracma::hessian(floglik, theta.til, data=data)
                 else -pracma::jacobian(fscore, theta.til, data=data)
	   score.til.data <- if(is.null(fscore)) pracma::grad(floglik,theta.til,data=data)
                         else fscore(theta.til,data)
       ## other main results from fit
       dpsi.dtheta <- pracma::grad(fpsi, theta.hat) 
       var.psi.hat <- dpsi.dtheta^2 * var.theta.hat  
       se.psi.hat <- sqrt(var.psi.hat)             
   	   ##### complete the observed informations
       j.hat.new <-  j.hat / dpsi.dtheta^2
   	  }
    r <-  sqrt(2 * (el.hat - el.til)) * sign(psi.hat - psival)
    if(ronly) out <- list(r=r, theta.hat=theta.hat, psi.hat = psi.hat, theta.hyp=theta.til, psi.hyp=psival)
    if(!ronly) 
     {
        ####### perform simulations 
    	if(trace) cat("start Monte Carlo computation","\n")
	    meanAll <- rep(0, 2) 
	    prodAll <- matrix(0, 2, 2) 
	    dataSim <- as.list(data,all.names=TRUE)
	    seed.in <- if(is.null(seed)) sample.int(2^30,1)
                   else seed  
	    set.seed(seed.in)
	    if(trace) pb <- txtProgressBar(style=3) 
	    for (i in 1:R)
         {          
           if( (i%%(R/10)==0) & trace) setTxtProgressBar(pb, i/R)
           dataSim <- datagen(theta.hat, data=data)  
           l1 <- floglik(theta.hat, dataSim)
           l0 <- floglik(theta.til, dataSim)
	         score.hat <-  if(is.null(fscore))  pracma::grad(f=floglik, x0=theta.hat, data=dataSim)    
	                     else fscore(theta.hat, dataSim)
	       ## get (psi,lam) scores from theta scores
           obj.score.new.hat = score.hat / dpsi.dtheta
           uhh <- c(obj.score.new.hat,  (l1-l0))    
           meanAll <- meanAll + uhh / R
           prodAll <- prodAll +  tcrossprod(uhh) / R             
          }
      if(trace) close(pb)  
      #######  computes SS approximations   
      covAll <- prodAll * R/(R-1) - tcrossprod(meanAll)  * R/(R-1)
      i.hat <- covAll[1, 1]       
      q <-  covAll[2, 1]        
      u <- q * sqrt(j.hat.new) / i.hat 
      #### final r* computation 
      INF <- (1/r) * log(u/r)   
      rs <- r + INF 
      out <- list(r=r,  INF=INF, rs=rs, theta.hat = theta.hat, se.theta.hat=se.theta.hat , 
        info.hat=j.hat, psi.hat = psi.hat, se.psi.hat = se.psi.hat,  
        theta.hyp=theta.til, psi.hyp=psival, seed=seed.in)
    }
  }  
  out$psidesc <- psidesc
  out$R <- R
  ## exit
  if( (!ronly) & (abs(r)<0.10))
    {
    	cat("Value under testing close to the MLE - there might be a singularity in r*\n")
    	warning("Value under testing close to the MLE - there might be a singularity in r*\n")
    }	
  return(structure(out, class="rstar"))
}


print.rstar <- function(x, digits = max(3, getOption("digits") - 3), ...) 
{ 
#
#  Print r* object
#
  cat("psi value under testing \n")
  print.default(format(x$psi.hyp, digits=digits), print.gap = 2, quote = FALSE)
  cat("Maximum likelihood estimate of psi \n")
  print.default(format(x$psi.hat, digits=digits), print.gap = 2, quote = FALSE)
  if(!is.null(x$se.psi.hat)) 
    {
     cat("Standard error of maximum likelihood estimate of psi \n")
     print.default(format(drop(x$se.psi.hat), digits=digits), print.gap = 2, quote = FALSE)
     }     
  cat("r statistic\n")
  print.default(format(x$r, digits=digits), print.gap = 2, quote = FALSE)
  if(!is.null(x$rs)) 
    {
     cat("r* statistic\n")
     print.default(format(drop(x$rs), digits=digits), print.gap = 2, quote = FALSE)
    }
  invisible(x)
}


summary.rstar <- function(object, ...){
#
#  Summary of r* testing
#
  digits <- max(3, getOption("digits") - 3)
  if(class(object)!="rstar")
    stop("\n'summary.rstar' designed for 'rstar' objects\n")
  if(!is.null(object$rs)) cat("\nTesting based on the r and r* statistics\n")
  else cat("\nTesting based on the r statistic\n")
  cat("-----------------------------------------------------------\n")
  if(!is.null(object$psidesc)) cat("Parameter of interest:        ",object$psidesc,"\n")
  else  cat("Parameter of interest:        User-defined function\n")
  if(!is.null(object$R))
    cat("Skovgaard covariances computed with", object$R,"Monte Carlo draws\n")  
  cat("psi value under testing: \n")
  print.default(format(object$psi.hyp, digits=digits), print.gap = 2, quote = FALSE)
  cat("-----------------------------------------------------------\n")
  cat("Estimates\n")
  if(!is.null(object$se.psi.hat)) 
    {
     cat("Maximum likelihood estimate of psi:\n")
     print.default(format(drop(object$psi.hat), digits=digits), print.gap = 2, quote = FALSE)
     cat("Standard error of maximum likelihood estimate of psi: \n")
     print.default(format(drop(object$se.psi.hat), digits=digits), print.gap = 2, quote = FALSE)
     }    
  else
  { 
    cat("Maximum likelihood estimate of psi:\n")
    out <- print.default(format(object$psi.hat, digits=digits), print.gap = 2, quote = FALSE)
  }
  cat("Maximum likelihood estimate of theta:\n")
  print.default(format(object$theta.hat, digits=digits), print.gap = 2, quote = FALSE)
  cat("Maximum likelihood estimate of theta under the null:\n")
  print.default(format(object$theta.hyp, digits=digits), print.gap = 2, quote = FALSE) 
  cat("-----------------------------------------------------------\n")
  cat("Test Statistics\n")
  if(!is.null(object$rs)) 
    {
     cat("Wald statistic P(r_wald<observed value; 1st order):\n")
     rw <- (object$psi.hat - object$psi.hyp) / drop(object$se.psi.hat)
     print.default(format(c(rw, pnorm(rw)), digits=digits), print.gap = 3, quote = FALSE)
    }
  cat("r statistic    P(r<observed value; 1st order):\n")
  print.default(format(c(object$r, pnorm(object$r)), digits=digits), print.gap = 3, quote = FALSE)
  if(!is.null(object$rs)) 
    {
     cat("r* statistic   P(r<observed value; 2nd order):\n")
     print.default(format(c(object$rs, pnorm(object$rs)), digits=digits), print.gap = 3, quote = FALSE)
     if(!is.null(object$NP))
       {
        cat("-----------------------------------------------------------\n")
        cat("Decomposition of high-order adjustment r*-r\n")
        cat("NP adjustment  INF adjustment:\n")
        print.default(format(c(object$NP, object$INF), digits=digits), print.gap = 3, quote = FALSE)
       }
    }    
  cat("-----------------------------------------------------------\n")
  invisible(object)
}




.rstar.stripped <- function(data, theta.hat, theta.til.in, psival, floglik, fpsi, datagen, fscore, R, seed, ronly,
                            j.hat, score.hat.data, constr.opt)
{
#
#  Internal function, r* without any embellishment
#  uses theta.hat rather than thetainit; seed is never NULL 
#
  p <- length(theta.hat)
  min.floglik <- function(theta, data) (-1) * floglik(theta,data)
  min.fscore <- if(is.null(fscore)) NULL else function(theta, data) (-1) * fscore(theta,data)
  el.hat <- floglik(theta.hat, data)
  if(length(theta.hat)>1)
   {
    psifcn.mod <- if (constr.opt=="solnp") function(theta, data) fpsi(theta)  	
                  else function(theta, data) fpsi(theta) - psival 
  	objHyp <- if (constr.opt=="solnp") solnp(theta.hat,fun=min.floglik, eqfun=psifcn.mod, eqB=psival,
	           control=list(trace=0), data=data)
	          else constrOptim.nl(theta.hat,fn=min.floglik, heq=psifcn.mod, gr=min.fscore,
	           control.outer=list(trace=FALSE), data=data) 
    theta.til <- objHyp$par    
    el.til <- floglik(theta.til,data)
   	psi.hat <- fpsi(theta.hat)
   	if(!ronly)
	   {
        j.til <-  if(is.null(fscore)) -pracma::hessian(floglik, theta.til, data=data)
                  else -pracma::jacobian(fscore, theta.til,data=data)
        score.til.data <- if(is.null(fscore)) pracma::grad(floglik,theta.til,data=data)
	                      else fscore(theta.til,data)
       }
	 r <-  sqrt(2 * (el.hat - el.til)) * sign(psi.hat - psival)
	 if(!ronly)
	   {
	    C.hat <- pracma::grad(fpsi, theta.hat)
	    C.til <- pracma::grad(fpsi, theta.til)
        k <- which(C.hat != 0)[1] 
        obj.info <- .newinfo(p, k, C.hat, C.til, j.hat, j.til, score.hat.data, score.til.data, theta.hat, theta.til, fpsi)
 	    j.hat.new <- obj.info$j.hat.new
	    j.til.new <- obj.info$j.til.new 
	   }
    if(ronly)  out<- list(r=r, psi.hyp=psival, theta.hyp=theta.til)
    if(!ronly) 
     {
        meanAll <- rep(0, 2*p+1) 
	    prodAll <- matrix(0, 2*p+1, 2*p+1)
	    dataSim <- as.list(data,all.names=TRUE)
	    set.seed(seed)
	    for (i in 1:R)
         { 
           dataSim <- datagen(theta.hat, data=data)  
           l1 <- floglik(theta.hat, dataSim)
           l0 <- floglik(theta.til, dataSim)
	       score.hat <-  if(is.null(fscore))  pracma::grad(f=floglik, x0=theta.hat, data=dataSim)    
	                     else fscore(theta.hat, dataSim)
	       score.til <-  if(is.null(fscore)) pracma::grad(f=floglik, x0=theta.til, data=dataSim) 
                          else fscore(theta.til, dataSim)
	       obj.score <- .newscores(p, k, C.hat, C.til, score.hat, score.til)
	       uhh <- c(obj.score$score.new.hat,  obj.score$score.new.til, (l1-l0))    
           meanAll <- meanAll + uhh / R
           prodAll <- prodAll +  tcrossprod(uhh) / R             
          }
      covAll <- prodAll * R/(R-1) - tcrossprod(meanAll)  * R / (R-1)
      S <-  covAll[1:p, (p+1):(2*p)] 
      i.hat <- covAll[1:p, 1:p]       
      i.hatInv <- qr.solve(i.hat,tol=10^-20)
      q <-  covAll[2*p+1, 1:p]        
      SS2 <- t(S) %*% i.hatInv %*% j.hat.new 
      SS1 <- q %*% i.hatInv %*% j.hat.new
      indpsi <- k  
      numU <- SS1[indpsi] -  SS2[-indpsi,indpsi] %*% solve(t(SS2[-indpsi,-indpsi])) %*% SS1[-indpsi] 
      j.hatInv.new <- solve(j.hat.new)  
      jProf.new <- 1 / j.hatInv.new[indpsi, indpsi] 
      u <- numU / sqrt(jProf.new)
      CPsi <- det(as.matrix(SS2[-indpsi, -indpsi])) / sqrt(det(as.matrix(j.til.new[-indpsi,-indpsi]   )) * det(as.matrix(j.hat.new[-indpsi, -indpsi])))  
      NP <- (1/r) * log(CPsi)  
      INF <- (1/r) * log(u/r)   
      rs <- r + NP + INF 
      out <- list(r=r, NP=NP, INF=INF, rs=rs, psi.hyp=psival, theta.hyp=theta.til)
     }
  } 
  if(length(theta.hat)==1)
   {
   	psifcn.mod <- function(theta, data) fpsi(theta) - psival   
   	objHyp <- nleqslv(theta.hat, psifcn.mod)
	theta.til <- objHyp$x
    el.til <- floglik(theta.til, data) 
    psi.hat <- fpsi(theta.hat)
    dpsi.dtheta <- pracma::grad(fpsi, theta.hat)    
    if(!ronly)
   	  { 
       j.til <-  if(is.null(fscore)) -pracma::hessian(floglik, theta.til, data=data)
                 else -pracma::jacobian(fscore, theta.til, data=data)
	   score.til.data <- if(is.null(fscore)) pracma::grad(floglik,theta.til,data=data)
                         else fscore(theta.til,data)
       }
	r <-  sqrt(2 * (el.hat - el.til)) * sign(psi.hat - psival)
	if(ronly) out <- list(r=r, psi.hyp=psival)
    if(!ronly) 
     {
     	meanAll <- rep(0, 2) 
	    prodAll <- matrix(0, 2, 2) 
	    dataSim <- as.list(data,all.names=TRUE)
	    set.seed(seed)
	    for (i in 1:R)
         { 
           dataSim <- datagen(theta.hat, data=data)  
           l1 <- floglik(theta.hat, dataSim)
           l0 <- floglik(theta.til, dataSim)
	       score.hat <-  if(is.null(fscore))  pracma::grad(f=floglik, x0=theta.hat, data=dataSim)    
	                     else fscore(theta.hat, dataSim)
	       obj.score.new.hat = score.hat / dpsi.dtheta
           uhh <- c(obj.score.new.hat,  (l1-l0))    
           meanAll <- meanAll + uhh / R
           prodAll <- prodAll +  tcrossprod(uhh) / R             
          }
      covAll <- prodAll * R/(R-1) - tcrossprod(meanAll)  * R/(R-1)
      i.hat <- covAll[1, 1]       
      q <-  covAll[2, 1]        
      j.hat.new <-  j.hat / dpsi.dtheta^2
      u <- q * sqrt(j.hat.new) / i.hat 
      INF <- (1/r) * log(u/r)   
      rs <- r + INF 
      out <- list(r=r,  INF=INF, rs=rs, psi.hyp=psival,  theta.hyp=theta.til)
    }
  }  
  ## exit
  return(out)
}


rstar.ci <- function(data, thetainit, floglik, fscore=NULL, fpsi, datagen, R=1000, seed=NULL,  
                     ronly=FALSE, psidesc=NULL,  constr.opt="solnp",
                     lower=NULL, upper=NULL, control=list(...), ...)
{
  
  ####### checking input arguments
  if(!is.list(data))
    {
     warning("data should be provided as a list\n")
     data <- as.list(data)
    } 
  if(!is.numeric(thetainit))
    stop("a starting point for the parameter theta is required \n")
  f0 <- floglik(thetainit, data)
  if(!is.numeric(f0))
     stop("problems in the loglikelihood function \n")
  if(!is.null(fscore))
   {
   	g0 <- fscore(thetainit, data)
    if(length(g0)!=length(thetainit))
      stop("size of starting point different from the size of score function \n")
   }
  if(!ronly)
    {
     data0 <- datagen(thetainit, data)
     f0 <- floglik(thetainit, data0)
     if(!is.numeric(f0))
     stop("problems in the function to simulate data \n")
    } 
   if( (constr.opt!="solnp") & (constr.opt!="alabama"))
   stop("constrained optimizer must be either 'solnp' or 'alabama'") 
    
   ###### full MLE, computed only once
   p <- length(thetainit)
   min.floglik <- function(theta, data) (-1) * floglik(theta,data)
   min.fscore <- if(is.null(fscore)) NULL else function(theta, data) (-1) * fscore(theta,data)
   obj.hat <- nlminb(thetainit, min.floglik, min.fscore, data=data)
   theta.hat <- obj.hat$par
   j.hat <- if(is.null(fscore)) -pracma::hessian(floglik, theta.hat, data=data)
           else  -pracma::jacobian(fscore, theta.hat, data=data)
   score.hat.data <- if(is.null(fscore)) pracma::grad(floglik, theta.hat, data=data)
	                else fscore(theta.hat, data=data)
   psi.hat <- fpsi(theta.hat)
   var.theta.hat <- solve(j.hat)
   dpsi.dtheta <- pracma::grad(fpsi, theta.hat) 
   var.psi.hat <- dpsi.dtheta %*% var.theta.hat %*% dpsi.dtheta
   se.psi.hat <- drop(sqrt(var.psi.hat)) 
   seed.in <- if(is.null(seed)) sample.int(2^30,1)
               else seed  
   ###### set up the computation            
   control <- do.call("rstar.ci.control", control)    
   away <- control$away
   npoints <- control$npoints
   trace <- control$trace
   stepsizefac <- control$stepsizefac             
   stepsize <- stepsizefac / npoints * se.psi.hat 
   maxstep <- control$maxstep
   theta.til.start <- theta.hat
   psivals.list <-  rvals.list <- NPvals.list <- INFvals.list <- c()
   psi.start <- psi.hat - se.psi.hat * away
   ###### if not requested from user, it gets the confidence limits by a stepwise approach
   if(is.null(lower) | is.null(upper))
   {
   	stepcount <- 0
   	repeat
   	  {
   	   if(trace) cat("Computing stats at",format(psi.start, digits=3),"\n")	
   	   rs.obj <- .rstar.stripped(data, theta.hat, theta.til.start, psi.start, floglik, fpsi, datagen, 
   	                             fscore, R, seed.in, ronly, j.hat, score.hat.data, constr.opt)  	
       psivals.list <- c(psi.start, psivals.list)
       rvals.list <- c(rs.obj$r, rvals.list)  
   	   val <- rs.obj$r
   	   if(!ronly)
         {
          if(length(thetainit)>1) NPvals.list <- c(rs.obj$NP, NPvals.list) 
          INFvals.list <- c(rs.obj$INF, INFvals.list)      	
          val <- rs.obj$rs
         }  
      theta.til.start <- rs.obj$theta.hyp    
   	  if(val>qnorm(0.995)) break 
   	  else  psi.start <- psi.start - stepsize
      stepcount <- stepcount + 1
   	  if(stepcount == maxstep) break  
   	 }
    lowerin <- psi.start 
    psi.start <- psi.hat + se.psi.hat * away
    stepcount <- 0
    repeat
   	 {
   	  if(trace) cat("Computing stats at",format(psi.start, digits=3),"\n") 	
   	  rs.obj <- .rstar.stripped(data, theta.hat, theta.til.start, psi.start, floglik, fpsi, datagen, 
   	                            fscore, R,  seed.in, ronly, j.hat, score.hat.data, constr.opt)  	
      psivals.list <- c(psivals.list, psi.start)
      rvals.list <- c(rvals.list, rs.obj$r)  
      val <- rs.obj$r
 	  if(!ronly)
          {
           if(length(thetainit)>1) NPvals.list <- c(NPvals.list, rs.obj$NP) 
           INFvals.list <- c(INFvals.list, rs.obj$INF)      	
           val <- rs.obj$rs
          }  
      theta.til.start <- rs.obj$theta.hyp    
   	  if(val<qnorm(0.005)) break 
   	  else  psi.start <- psi.start + stepsize
   	  stepcount <- stepcount + 1
   	  if(stepcount == maxstep) break
   	 }  
    upperin <- psi.start
    psi.out <- psivals.list
    r.out <- rvals.list
    if((!ronly) & (length(thetainit)==1))   INF.out <-  INFvals.list
    if((!ronly) & (length(thetainit)>1))   {INF.out <-  INFvals.list; NP.out <- NPvals.list}
   } 
   #########  only the lower-upper case, grid of psi values is available
   if( (!is.null(lower)) &  (!is.null(upper)) )
   {
    psi.left.out <-  lower
    psi.left.in <-  psi.hat - away * se.psi.hat 
    psi.right.in <- psi.hat + away * se.psi.hat  
    psi.right.out <- upper
    if(psi.left.in >= psi.left.out)
      warning("lower value too high\n")
    if(psi.right.in >= psi.right.out)
      warning("upper value too low\n")
    psi.grid.left <- seq(psi.left.in, psi.left.out, l=npoints)  
    psi.grid.right <- seq(psi.right.in, psi.right.out, l=npoints)  
    if(length(thetainit)>1)  NP.out <- rep(NA, 2*npoints)
    INF.out <- r.out <- psi.out <- rep(NA, 2*npoints)
    theta.til.start <- theta.hat
    for(i in 1:length(psi.grid.left))
      {
       psival <- psi.grid.left[i]
       if(trace) cat("Computing stats at",format(psival, digits=3),"\n") 
       rs.obj <- .rstar.stripped(data, theta.hat, theta.til.start, psival, floglik, fpsi, datagen, fscore, 
       	                          R, seed.in, ronly, j.hat, score.hat.data,  constr.opt)
       r.out[i] <- rs.obj$r      
       psi.out[i] <- rs.obj$psi.hyp                     
       if(!ronly)
         {
          if(length(thetainit)>1) NP.out[i] <- rs.obj$NP
          INF.out[i] <- rs.obj$INF	
         }  
      theta.til.start <- rs.obj$theta.hyp
    }
  theta.til.start <- theta.hat  
  for(i in 1:length(psi.grid.right))
     {
      psival <- psi.grid.right[i]
      if(trace) cat("Computing stats at",format(psival, digits=3),"\n")
      rs.obj <- .rstar.stripped(data, theta.hat, theta.til.start, psival, floglik, fpsi, datagen, fscore, 
       	                          R, seed.in, ronly, j.hat, score.hat.data,  constr.opt)
       	                          
      r.out[i+npoints] <- rs.obj$r      
      psi.out[i+npoints] <- rs.obj$psi.hyp                     
      if(!ronly)
        {
         if(length(thetainit)>1) NP.out[i+npoints] <- rs.obj$NP
         INF.out[i+npoints] <- rs.obj$INF	
        }  
      theta.til.start <- rs.obj$theta.hyp
     }
   }
   ###### lines up the results
   ord <- order(psi.out)
   if(ronly) out <- list(psivals = psi.out[ord], rvals=r.out[ord])
   if((!ronly) & (length(thetainit)==1)) out <- list(psivals = psi.out[ord], rvals=r.out[ord], 
                                  INFvals=INF.out[ord], rsvals=INF.out[ord]+r.out[ord])
   if((!ronly) & (length(thetainit)>1)) out <- list(psivals = psi.out[ord], rvals=r.out[ord], 
                                  INFvals=INF.out[ord], NPvals=NP.out[ord], rsvals=INF.out[ord]+r.out[ord]+NP.out[ord])       
   ###### gets confidence intervals
   obj.r <- smooth.spline(out$rvals, out$psivals, all.knots=TRUE)
   if(!ronly) obj.rs <- smooth.spline(out$rsvals, out$psivals,  all.knots=TRUE)
   alpha <-  c(0.10, 0.05, 0.01)
   CIr <- matrix(0, nrow=3, ncol=2)
   if(!ronly) CIrs <- matrix(0, nrow=3, ncol=2)
   for (i in 1:3)
     {
      CIr[i,] <- sort(predict(obj.r, c(qnorm(alpha[i]/2), qnorm(1-alpha[i]/2)))$y) 
      if(!ronly) CIrs[i,] <- sort(predict(obj.rs, c(qnorm(alpha[i]/2), qnorm(1-alpha[i]/2)))$y) 
     }
   out <- c(out, list(CIr=CIr))
   if(!ronly)  out <- c(out, list(CIrs=CIrs))
   out$psidesc <- psidesc
   out$R <- R
   out$seed <- seed.in
   ###### exit
   return(structure(out, class="rstarci"))    
}


rstar.ci.control <- function(npoints=10, away=0.3, stepsizefac=3, maxstep=50, trace=TRUE) 
{
  list(away=away, npoints=npoints, stepsizefac=stepsizefac, maxstep=maxstep, trace=trace) 
}  


print.rstarci <- function(x, digits = max(3, getOption("digits") - 3), ...){
#
#  Print confidence intervals
#
    cat("Confidence interval calculations based on likelihood asymptotics\n")
    cat("1st-order\n")
    cat("         90%                         95%                         99%     \n")
    rint90 <-  format(x$CIr[1,], digits=digits)
    rint95 <-  format(x$CIr[2,], digits=digits)
    rint99 <-  format(x$CIr[3,], digits=digits)
    cat("(",rint90[1]," , ",rint90[2],")         (",rint95[1]," , ",rint95[2],")         (",rint99[1]," , ",rint99[2],")\n")
    if(!is.null(x$CIrs))
       {
       	cat("2nd-order\n")
        cat("         90%                         95%                         99%     \n")
        rsint90 <-  format(x$CIrs[1,], digits=digits)
        rsint95 <-  format(x$CIrs[2,], digits=digits)
        rsint99 <-  format(x$CIrs[3,], digits=digits)
        cat("(",rsint90[1]," , ",rsint90[2],")        (",rsint95[1]," , ",rsint95[2],")        (",rsint99[1]," , ",rsint99[2],")\n")
       }
    invisible(x)
}


summary.rstarci <- function(object, ...){
#
#  Summary of confidence intervals
#
    if(class(object)!="rstarci")
       stop("\n'summary.rstarci' designed for 'rstarci' objects\n")
    digits <- max(3, getOption("digits") - 3)
    cat("Confidence interval calculations based on likelihood asymptotics\n")
    cat("-----------------------------------------------------------------------------\n")
    if(!is.null(object$psidesc)) cat("Parameter of interest:        ",object$psidesc,"\n")
    else  cat("Parameter of interest:        User-defined function\n")
    cat("Calculations based on a grid of", length(object$psivals),"points\n") 
    if(!is.null(object$R))
    cat("Skovgaard covariances computed with", object$R,"Monte Carlo draws\n")  
    cat("-----------------------------------------------------------------------------\n")
    cat("1st-order\n")
    cat("         90%                         95%                         99%     \n")
    rint90 <-  format(object$CIr[1,], digits=digits)
    rint95 <-  format(object$CIr[2,], digits=digits)
    rint99 <-  format(object$CIr[3,], digits=digits)
    cat("(",rint90[1]," , ",rint90[2],")       (",rint95[1]," , ",rint95[2],")       (",rint99[1]," , ",rint99[2],")\n")
    if(!is.null(object$CIrs))
       {
       	cat("2nd-order\n")
        cat("         90%                         95%                         99%     \n")
        rsint90 <-  format(object$CIrs[1,], digits=digits)
        rsint95 <-  format(object$CIrs[2,], digits=digits)
        rsint99 <-  format(object$CIrs[3,], digits=digits)
        cat("(",rsint90[1]," , ",rsint90[2],")       (",rsint95[1]," , ",rsint95[2],")       (",rsint99[1]," , ",rsint99[2],")\n")
        if(!is.null(object$NPvals))
          {
           cat("-----------------------------------------------------------------------------\n")
           cat("Decomposition of high-order adjustment\n")
           cat("Nuisance parameter adjustment (NP)\n")
           sNP <- format(summary(object$NPvals,digits=digits))
           print.default(sNP, print.gap = 2, quote = FALSE)	
           cat("Information adjustment (INF)\n")
           sINF <- format(summary(object$INFvals,digits=digits))
           print.default(sINF, print.gap = 2, quote = FALSE)
         }
       }      
   cat("-----------------------------------------------------------------------------\n")    
   invisible(object)
}


plot.rstarci <- function(x,colrs=2,ltyrs=1,...)
{
  if(class(x)!="rstarci") stop("\n'plot.rstarci' designed for 'rstarci' objects\n")
  yrange <-  if(!is.null(x$rsvals)) c( min(c(-3, x$rvals, x$rsvals)), max( c(3, x$rvals, x$rsvals)))
             else   c(min(c(-3, x$rvals)), max(c(3, x$rvals)))
  mar <- c(5, 4, 4, 5) + 0.1
  par(mar=mar)
  xlab <- if(is.null(x$psidesc)) expression(psi) else x$psidesc
  plot(x$psivals, x$rvals, type="n", ylab="r values", ylim=yrange, xlab=xlab, las=1)
  grid(NA, 20, lwd = 1.5)
  lines(spline(x$psivals, x$rvals))
  if(!is.null(x$rsvals)) lines(spline(x$psivals, x$rsvals), col=colrs, lty=ltyrs)
  points(x$psivals, x$rvals, pch=16)
  if(!is.null(x$rsvals)) points(x$psivals, x$rsvals, col=colrs, pch=16)
  rlab <- c(0.001, 0.01, 0.05, 0.10 ,0.5, 0.10, 0.05, 0.01, 0.001)
  rpos <- qnorm(c(rlab[1:5], 1-rlab[6:9]))
  axis(4, at=rpos, labels=c("0.001", "0.01", "0.05", "0.10" ,"0.5", "0.10", "0.05", "0.01", "0.001"), las=1)
  mtext("P values",side=4,line=3)
  if(!is.null(x$rsvals)) legend("topright",col=c(1,colrs),lty=c(1,ltyrs), lwd=1.5, legend=c("1st order","2nd order"), bty="n", cex=1.25)
  else  legend("topright",col=c(1),lty=1, lwd=1.5, legend=c("1st order"), bty="n", cex=1.25)       
  invisible(x)
}


.mpl <- function(data, mle, psival, floglik, fscore, indpsi, datagen, R, seed, plonly, trace)
{
#
#  mpl computation; internal function
#  
  ####### checking input arguments
  if(!is.list(data))
    {
     warning("data should be provided as a list\n")
     data <- as.list(data)
    } 
  if ((!is.numeric(mle)) & (plonly) )
    stop("initial value of the parameter theta is required \n")
  if ((!is.numeric(mle)) & (!plonly) )
    stop("mle of the parameter theta is required \n")
  if(!is.numeric(psival))
    stop("a value for the parameter of interest is required \n")
  f0 <- floglik(mle, data)
  if(!is.numeric(f0))
     stop("problems in the loglikelihood function \n")
  if(!is.null(fscore))
   {
   	g0 <- fscore(mle, data)
   	if(length(g0)!=length(mle))
      stop("size of starting point different from the size of score function \n")
   }
  if(!plonly)
    {
     data0 <- datagen(mle, data)
     f0 <- floglik(mle, data0)
     if(!is.numeric(f0))
     stop("problems in the function to simulate data \n")
    }
  if(!is.numeric(indpsi))
     stop(paste("index for psi must be a set of numbers in the range 1-",length(mle)))
  if( (min(indpsi)<1)| (max(indpsi)>length(mle)) ) 
     stop(paste("index for psi must be a set of numbers in the range 1-",length(mle))) 
  
  
  ####### getting the constrained MLE
  if(trace) cat("Computing stats at",format(psival, digits=5),"\n")
  p <- length(mle)
  theta.hat <- mle
  if(trace) cat("get mle under the null....","\n")
  
  min.floglik <-function(lambda)
	{ 
	  theta <- rep(0,p)
	  theta[indpsi] <- psival
      theta[-indpsi] <- lambda
	  out <- floglik(theta,data) * (-1)
	  return(out)
	}
	 
  min.fscore <- if(is.null(fscore)) NULL
  else function(lambda)
	       { 
	       	 theta<- rep(0,p)
	         theta[indpsi]<- psival
             theta[-indpsi]<- lambda
	         out<- fscore(theta,data)[-indpsi] * (-1)
	  return(out)
	} 
  init <- mle[-indpsi]
  objHyp <- nlminb(as.numeric(init), min.floglik, min.fscore)
  theta.til <- theta.hat
  theta.til[indpsi] <- psival
  theta.til[!indpsi] <- objHyp$par
  el.til <- floglik(theta.til, data)
  if(plonly) out <- el.til
  if(!plonly) 
     {
        j.til <-  if(is.null(fscore)) -pracma::hessian(floglik, theta.til, data=data)                  
                  else -pracma::jacobian(fscore, theta.til, data=data)
        j.hat <-  if(is.null(fscore)) -pracma::hessian(floglik, theta.hat, data=data)                  
                  else -pracma::jacobian(fscore, theta.hat, data=data)          
        ####### perform simulations 
    	if(trace) cat("start Monte Carlo computation","\n")
	    meanAll <- rep(0, 2*p) 
	    prodAll <- matrix(0, 2*p, 2*p)
	    dataSim <- as.list(data,all.names=TRUE)
	    if(!is.null(seed)) set.seed(seed)
	    if(trace) pb <- txtProgressBar(style=3) 
	    for (i in 1:R)
         { 
           if( (i%%(R/10)==0) & trace) setTxtProgressBar(pb, i/R)
           dataSim <- datagen(theta.hat, data=data)  
           score.hat <-  if(is.null(fscore))  pracma::grad(f=floglik, x0=theta.hat, data=dataSim)    
	                     else fscore(theta.hat, dataSim)
	       score.til <-  if(is.null(fscore)) pracma::grad(f=floglik, x0=theta.til, data=dataSim) 
                          else fscore(theta.til, dataSim)
	       uhh <- c(score.hat, score.til)    
           meanAll <- meanAll + uhh / R
           prodAll <- prodAll +  tcrossprod(uhh) / R             
        }
      if(trace) close(pb)   
      #######  computes SS approximations 
      covAll <- prodAll * R/(R-1) - tcrossprod(meanAll)  * R / (R-1)
      S <-  covAll[1:p, (p+1):(2*p)] 
      i.hat <- covAll[1:p, 1:p]       
      i.hatInv <- qr.solve(i.hat, tol=10^-20)
      SS2 <- t(S) %*% i.hatInv %*% j.hat
      out <-  el.til + 0.5 * log(det(j.til[-indpsi, -indpsi])) - log(det(SS2[-indpsi, -indpsi]))
  }
  ## exit
  if(trace) cat("Function value",format(out, digits=5),"\n\n")
  return(out)  
}


logPL <- function(psival, data, thetainit, floglik, fscore=NULL, indpsi, minus=FALSE, trace=FALSE)
{
#	
# Wrapper function that computes minus the profile likelihood 
# R and seed are not used
#
  out <-  .mpl(data=data, mle=thetainit, psival=psival, floglik=floglik, indpsi=indpsi, datagen=NULL, fscore=fscore, R=NULL, 
  seed=NULL, plonly=TRUE, trace=trace)
  flag <- as.numeric(!minus) - as.numeric(minus) 
  return(out * flag)
}	


logMPL <- function(psival, data, mle, floglik, fscore=NULL, indpsi, datagen, R=500, seed=NULL, minus=FALSE, trace=FALSE)
{
#	
# Wrapper function that computes minus the modified profile likelihood
#
  out <- .mpl(data=data, mle=mle, psival=psival, floglik=floglik, fscore=fscore, indpsi=indpsi, datagen=datagen,
              R=R, seed=seed, trace=trace, plonly=FALSE)
  flag <- as.numeric(!minus) - as.numeric(minus)             
  return(out * flag)
}	
	

	







