## file nlreg/R/nlreg.R, v 1.2-2 2014-03-31
##
##  Copyright (C) 2000-2014 Ruggero Bellio & Alessandra R. Brazzale 
##
##  This file is part of the "nlreg" package for R.  This program is 
##  free software; you can redistribute it and/or modify it under the 
##  terms of the GNU General Public License as published by the Free 
##  Software Foundation; either version 2 of the License, or (at your 
##  option) any later version.
##
##  This library is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
##  02111-1307 USA or look up the web page
##  http://www.gnu.org/copyleft/gpl.html.
##
##  Please send any comments, suggestions or errors found to:
##  Alessandra R. Brazzale, Department of Statistics, University of
##  Padova, Via C. Battisti 241/243, 35121 Padova (PD), Italy.
##  Email: alessandra.brazzale@unipd.it
##  Web: http://www.stat.unipd.it/~brazzale

Dmean <- function(nlregObj, hessian = TRUE)       
{
  rc <- nlregObj$coef 
  do.call("deriv3", list(expr = nlregObj$meanFun, namevec = names(rc),
	                 function.arg = c(names(nlregObj$ws$allPar), "..."), 
                         hessian = hessian))
}

Dvar <- function(nlregObj, hessian = TRUE)
{
  rc <- nlregObj$coef 
  vp <- nlregObj$varPar
  namevec <- if(nlregObj$ws$xVar)
               c(names(rc), names(vp))
             else names(vp)
  do.call("deriv3", list(expr = nlregObj$varFun, namevec = namevec,
                         function.arg = c(names(nlregObj$ws$allPar), "..."), 
                         hessian = hessian))
}  

nlreg <- function(formula, weights = NULL, 
                  data = sys.frame(sys.parent()), start, 
                  offset = NULL, subset = NULL, 
                  control = list(x.tol = 1e-6, rel.tol = 1e-6,
                                 step.min = 1/2048, maxit = 100),
                  trace = FALSE, hoa = FALSE)
{
  m <- match.call()  
  if( missing(formula) )  
    stop("formula is missing, with no default")
  if( !is.null(weights) & !is.call(weights) )
    stop("variance function must be a call")
  if( missing(start) )  
    stop("starting values are missing, with no default")
  if( !is.null(offset) & is.null(names(offset)) )
    stop("\'offset\' argument must be a named vector")
  mFormula <- formula
  mFormula[[2]] <- NULL
  allx <- c(all.vars(formula), if(is.call(weights)) all.vars(weights))
  if( !all(match(names(start), allx, nomatch=FALSE)) )
    stop("some starting values are missing, with no default")
  if( !is.null(offset) )
  {
    if( (names(offset) != "logs") && 
        (!match(names(offset), allx, nomatch=FALSE)) )
      stop("offset variable is not included in original formula")
    if( match(names(offset), names(start), nomatch=FALSE) )
      stop("no starting value admitted for offset parameter")
  }
  missingData <- missing(data)
  if( !missingData )
    if( !is.data.frame(data) ) 
      stop("\"data\" is not a data frame")
  if( !missing(control) )
  {
    if( is.null(control$x.tol) )
      control$x.tol <-  1e-6
    if( is.null(control$rel.tol) )
      control$rel.tol <-  1e-6
    if( is.null(control$step.min) )
      control$step.min <-  1/2048
    if( is.null(control$maxit) )
      control$maxit <- 100
  }	
  if( is.null(weights) )
  {
    if(trace)
      cat("\nconstant variance function: switch to \"nls\"\n\n")
    if( !missingData )
    {
     new.data <- data
      if( !is.null(offset) && (names(offset) != "logs") )
      {
        new.data <- cbind(new.data, offset)
        names(new.data)[dim(new.data)[[2]]] <- names(offset)
      }
      fit.temp <- list(as.name("nls"), 
                       formula=formula, data=as.name("new.data"), 
                       start=as.name("start"),
                       control=nls.control(maxiter=control$maxit,
                                          tol=control$rel.tol,
                                          minFactor=control$step.min),
                       algorithm="default", trace=trace)
      if( !is.null(subset) )
        fit.temp <- c(fit.temp, subset=subset)
      fit.temp <- eval(as.call(fit.temp))
    }
    else
    {
      if( !is.null(offset) )
        do.call("<-", list(names(offset), offset))
      fit.temp <- list(as.name("nls"), 
                       formula=formula, start=as.name("start"),
                       control=nls.control(maxiter=control$maxit,
                                           tol=control$rel.tol,
                                           minFactor=control$step.min),
                       algorithm="default", trace=trace)
      if( !is.null(subset) )
        fit.temp <- c(fit.temp, subset=subset)
      fit.temp <- eval(as.call(fit.temp))
    }
    n <- length(residuals(fit.temp))
    logs <- if( is.null(offset) || (names(offset) != "logs") )
              log(sum(residuals(fit.temp)^2)/n)
            else offset
    if( !missingData )
    {
      .nd <- names(data)
      .md <- match(".addVar", .nd, nomatch=0)
      if(.md > 0)
      {
        .nd[.md] <- NA	        
        cov.pos <- as.logical(match(.nd, all.vars(formula[[3]]), nomatch=0))
        xx <- unlist(data[cov.pos])
        yval <- eval(formula[[2]], data)
        new.data <- list( xx=xx, repl=rep(1, n), dupl=1:n, t1=yval, t2=yval^2)
        new.data <- lapply(new.data, as.vector)
        names(new.data) <- c(names(data)[cov.pos], "repl", "dupl", "t1", "t2")
      }
      else
      {
        cov.pos <- match(names(data), all.vars(mFormula))
        cov.pos <- !is.na(cov.pos)    
        new.data <- list( (xx <- unique(unlist(data[cov.pos]))),  
                  tapply(unlist(data[cov.pos]), unlist(data[cov.pos]),
                         function(x) sum(duplicated(x)) + 1)[order(xx)],
                  match(unlist(data[cov.pos]), 
                        unique(unlist(data[cov.pos]))),
                  tapply(eval(formula[[2]], data), 
                         unlist(data[cov.pos]),
                         function(x) sum(x))[order(xx)],
                  tapply(eval(formula[[2]], data), 
                         unlist(data[cov.pos]),
                         function(x) sum(x^2))[order(xx)] )          
        new.data <- lapply(new.data, as.vector)
        names(new.data) <- c(names(data)[cov.pos], "repl", "dupl", 
                             "t1", "t2")
      }                       
    }
    else
    {
      yval <- eval(formula[[2]])
      new.data <- list( repl=rep(1, n), dupl=1:n, t1=yval, t2=yval^2)
      new.data <- lapply(new.data, as.vector)
      names(new.data) <- c("repl", "dupl", "t1", "t2")
    }
    fit <- list( coef      = coef(fit.temp),
                 varPar    = if( is.null(offset) || 
                                 (names(offset) != "logs") ) 
                               c(logs=logs) else NULL,
                 offset    = offset,	  
                 logLik    = -n/2*logs - n/2*log(2*pi) - 
                               sum(residuals(fit.temp)^2)/exp(logs)/2,
                 meanFun   = mFormula,
                 varFun    = as.formula(call("~", 
                                       call("exp", as.name("logs")))),
                 data      = new.data,
                 fitted    = fitted(fit.temp),
                 weights   = rep(exp(logs), n),
                 residuals = residuals(fit.temp)/sqrt(exp(logs)), 
                 start     = start,
                 call      = m )
    allpar <- if( is.null(offset) )
                c(fit$coef, fit$varPar)
              else c(fit$coef, fit$varPar, offset)
    fit <- c(fit, list(ws = list( allPar      = allpar, 
                                  homVar      = TRUE, 
                                  xVar        = FALSE,
                                  hoa         = hoa,
                                  missingData = missingData,	
                                  frame       = if(!missingData) 
                                                 m$data else NULL ) ))
    attr(fit, "class") <- c("nlreg", "nls")
    if(hoa)
    {
      if( !any(match("nlreg.diag", unlist(lapply(sys.calls(), 
                       function(x) x[[1]])), nomatch=FALSE)) )
        cat("\ndifferentiating mean and variance function -- may take a while\n")
      fit$ws$md <- Dmean(fit)
      fit$ws$vd <- if( is.null(offset) || (names(offset) != "logs") ) 
                     Dvar(fit) else NULL
    }
  }
  else
  {			
    minRegCoef <- function(regCoef, varPar, offset, label, meanFun, 
                           varFun, data, subset, isCoef, lastIter)
    {
      new.data <- data
      if( is.null(new.data) )
      { 
        new.frame <- as.list(regCoef)
        names(new.frame) <- if(isCoef) label else names(regCoef)
      } 
      else
      {
        class(new.data) <- "list"
        tmp <- regCoef
        names(tmp) <- if(isCoef) label else names(regCoef)
        new.frame <- c(new.data, tmp)
      }
      if( !is.null(varPar) )
      {
        tmp <- varPar
        names(tmp) <- if(isCoef) names(varPar) else label
        new.frame <- c(new.frame, tmp) 
      }
      if( !is.null(offset) )
      {
        tmp <- offset
        names(tmp) <- names(offset)
        new.frame <- c(new.frame, tmp) 
      }
      num <- eval( call("-", meanFun[[2]], meanFun[[3]]), 
                   envir=new.frame )
      num <- num^2
      den <- eval( varFun[[2]], envir=new.frame )
      if( !is.null(subset) )
      {
        num <- num[subset] ; den <- den[subset]
      }
      nulloff <- is.null(offset)
      if( !nulloff )  nulloff <- (names(offset) != "logs")
      res <- if( nulloff )
             {
               if(!lastIter) 
                 log(sum(num/den)/length(num))*length(num) + 
                   sum(log(den)) + length(num) + 
                   length(num)*log(2*pi)
               else  log(sum(num/den)/length(num))
             }
             else
             {
               sigma2 <- exp(offset)
               if(!lastIter)
                 sum(log(den)) + sum(num/den)/sigma2 + 
                   length(num)*log(sigma2) + length(num)*log(2*pi)
               else offset
             }            
      res
    }
    minVarPar <- function(varPar, regCoef, offset, label, meanFun, 
                          varFun, data, subset)
    {
      res <- minRegCoef(regCoef=regCoef, varPar=varPar, offset=offset,
                        label=label, meanFun=meanFun, varFun=varFun, 
                        data=data, subset=subset, isCoef=FALSE, 
                        lastIter=FALSE)
      res
    }        
    fitted.val <- function(formula, weights, par, data, subset)
    {
      new.data <- data
      if( is.null(data) )
      { 
        new.frame <- as.list(par)
        names(new.frame) <- names(par)
      } 
      else
      {
        class(new.data) <- "list"
        tmp <- as.list(par)
        names(tmp) <- names(par)
        new.frame <- c(new.data, tmp)
      }
      fitted <- eval(formula[[3]], envir=new.frame)
      weights <- eval(weights[[2]], envir=new.frame)
      resid <- eval( call("-", formula[[2]], formula[[3]]), 
                     envir=new.frame )/sqrt(weights)
      if( !is.null(subset) )
      {
        fitted <- fitted[subset]
        weigths <- weights[subset]
        resid <- resid[subset]
      }
      list(fitted = fitted, weights = weights, resid=resid)
    }
    para <- names(start)
    coef.idx <- match(para, all.vars(formula[[3]]), nomatch=0)
    var.idx <- (coef.idx == 0)        
    coef.idx <- (coef.idx != 0)      
    tol <- control$x.tol ; test1 <- 0 ; test2 <- 0
    reltol <- control$rel.tol
    maxit <- control$maxit
    go <- TRUE
    coef.old <- coef.new <- if( any(coef.idx) )    
                              unlist(start[coef.idx]) else NULL
    var.old <- var.new <- if( any(var.idx) ) 
                            unlist(start[var.idx]) else NULL
    idx <- 0
    if(trace)  cat("\n")
    while(go) 
    {
      if( any(coef.idx) )
      {                  
	coef.new <- optim(par=coef.old, fn=minRegCoef, varPar=var.old,
                          offset=offset, label=names(coef.old), 
                          meanFun=formula, varFun=weights, 
                          subset=subset,
                          data=if(!missingData) data else NULL, 
                          isCoef=TRUE, lastIter=FALSE, 
#                          control=control, method="BFGS")
                          control=list(maxit=control$maxit, reltol=control$rel.tol), 
                          method="BFGS")
        log.lik <- -coef.new$value/2
        coef.new <- coef.new$par
        test2 <- max(abs((coef.old-coef.new)/coef.new)) 
        coef.old <- coef.new
      }            
      if( any(var.idx) )
      {
        var.new <- optim(par=var.old, fn=minVarPar, regCoef=coef.old, 
                         offset=offset, label=names(var.old), 
                         meanFun=formula, varFun=weights, 
                         subset=subset, 
                         data=if(!missingData) data else NULL, 
#                         control=control, method="BFGS")
                         control=list(maxit=control$maxit, reltol=control$rel.tol), 
                         method="BFGS") 
        log.lik <- -var.new$value/2 
        var.new <- var.new$par
        test1 <- max(abs((var.old-var.new)/var.new)) 
        var.old <- var.new
      }
      if( idx==0 )
      {
        new.lik <- log.lik
        test3 <- reltol + 1
      }
      else 
      {
        old.lik <- new.lik
        new.lik <- log.lik
        test3 <- max(abs((old.lik-new.lik)/new.lik))
      }
      idx <- idx + 1
      if(trace)
        cat(paste("iteration", idx, ": log likelihood =", 
                  format(log.lik), "\n"))
      if ( ((test1 < tol) &&  (test2 < tol)) || (test3 < reltol) || 
           (idx == maxit) )
      go <- FALSE
    }
    if( (maxit > 1) && (idx == maxit) )
      warning(paste("\nconvergence not obtained in", maxit, 
                    "iterations"))
    coeff <- if( any(coef.idx) ) coef.old else NULL
    nulloff <- is.null(offset)
    if(!nulloff)  nulloff <- (names(offset) != "logs")
    varPar <- if( nulloff )
                  c(logs=minRegCoef(coef.old, var.old, offset=offset,
                                    names(coef.old), formula, weights,
                                    data=if(!missing(data)) 
                                           data else NULL, 
                                    subset=subset, TRUE, TRUE))
              else NULL
    if( any(var.idx) )  varPar <- c(var.old, varPar) 
    allpar <- if( is.null(offset) )
                c(coeff, varPar)
              else c(coeff, varPar, offset)
    varFun <- as.formula(call("~", call("*", weights[[2]], 
                              call("exp", as.name("logs")))))
    more.fit <- fitted.val(formula, varFun, allpar, 
                           data=if(!missingData) data else NULL, 
                           subset=subset)	
    n <- length(more.fit$resid)
    if( !missing(data) )
    {    	
      .nd <- names(data)
      .md <- match(".addVar", .nd, nomatch=0)
      if(.md > 0)
      {
        .nd[.md] <- NA	        
        cov.pos <- as.logical(match(.nd, all.vars(formula[[3]]), nomatch=0))
        xx <- unlist(data[cov.pos])
        yval <- eval(formula[[2]], data)
        new.data <- list( xx=xx, repl=rep(1, n), 1:n, t1=yval, t2=yval^2)
        new.data <- lapply(new.data, as.vector)
        names(new.data) <- c(names(data)[cov.pos], "repl", "dupl", "t1", "t2")
      }
      else
      {
        cov.pos <- match(names(data), all.vars(formula[[3]]))		## 20.05.13      
        cov.pos <- !is.na(cov.pos) 
        new.data <- list( (xx <- unique(unlist(data[cov.pos]))),  
             tapply(unlist(data[cov.pos]), unlist(data[cov.pos]),
                    function(x) sum(duplicated(x)) + 1)[order(xx)], 
             match(unlist(data[cov.pos]), 
                   unique(unlist(data[cov.pos]))),	
             tapply(eval(formula[[2]], data), unlist(data[cov.pos]),
                    function(x) sum(x))[order(xx)],
             tapply(eval(formula[[2]], data), unlist(data[cov.pos]),
                    function(x) sum(x^2))[order(xx)] )    
        new.data <- lapply(new.data, as.vector)
        names(new.data) <- c(names(data)[cov.pos], "repl", "dupl", 
                             "t1", "t2")
      }
    }
    else
    {
      yval <- eval(formula[[2]])
      new.data <- list( repl=rep(1, n), 1:n, t1=yval, t2=yval^2)
      new.data <- lapply(new.data, as.vector)
      names(new.data) <- c("repl", "dupl", "t1", "t2")
    }
    xVar <- any( match(names(coeff), all.vars(weights[[2]]), 
                       nomatch=0) )
    fit <- list( coef      = coeff,  
                 varPar    = varPar,
                 offset    = offset,
                 logLik    = - minRegCoef(coef.old, var.old, 
                                          offset=offset, 
                                          names(coef.old), formula, 
                                          weights, 
                                          data=if(!missing(data)) 
                                                 data 
                                               else NULL, 
                                          subset, TRUE, FALSE)/2,
                 meanFun   = as.formula(call("~", formula[[3]])),
                 varFun    = varFun,
                 data      = new.data, 
                 fitted    = more.fit$fitted,
                 weights   = more.fit$weights,
                 residuals = more.fit$resid, 
                 start     = start,
                 call      = m,
                 ws        = list( allPar      = allpar, 
                                   homVar      = FALSE, 
                                   xVar        = xVar,
                                   hoa         = hoa,
                                   missingData = missingData,
                                   frame       = if(!missingData) 
                                                   m$data 
                                                 else NULL,
                                   iter        = idx ) )	
    attr(fit, "class") <- c("nlreg", "nls")
    if(hoa)
    {
      if( !any(match("nlreg.diag", unlist(lapply(sys.calls(), 
                     function(x) x[[1]])), nomatch=FALSE)) )
      cat("\ndifferentiating mean and variance function -- may take a while\n")
      fit$ws$md <- Dmean(fit)
      fit$ws$vd <- Dvar(fit)		
    }
  }
  fit
}

logLik.nlreg <- function(object, ...)
{
  if( length(list(...)) )
    warning("extra arguments discarded")
  nobs <- length(object$residuals)      
  coef <- object$coef
  varPar <- object$varPar
  p <- length(coef) + length(varPar)
  val <- object$logLik
  attr(val, "nobs") <- nobs
  attr(val, "npar") <- p 
  attr(val, "df") <- nobs - p 
  class(val) <- "logLik"
  val
}

print.nlreg <- function(x, digits = max(3, getOption("digits")-3), 
                        ...)
{
  if(!is.null(cl <- x$call)) 
  {
    cat("Formula:\n")
    dput(x$call$formula)
    cat("\nVariance function:\n")
    dput(x$varFun)		
  }
  coef <- x$coef
  if( !is.null(coef) )
  {
    cat("\nRegression coefficients:\n")
    print(coef, digits=digits, ...)
  }
  else cat("\nNo regression coefficient\n")
  varPar <- x$varPar
  if( !is.null(varPar) )
  {
    cat("\nVariance parameters:\n")
    print(varPar, digits=digits, ...)
  }
  else cat("\nNo variance parameter\n")
  of <- x$offset
  if( !is.null(of) )
  {
    cat("\nParameter of interest:\n")
    print(of, digits=digits, ...)
  }
  else cat("\nNo interest parameter\n")	
  nobs <- length(x$residuals)      
  cat("\nTotal number of observations:", nobs)
  cat("\nTotal number of parameters:", 
      sum(length(coef)+length(varPar)))
  cat("\n-2*Log Likelihood", 
      format(-2 * x$logLik, digits=digits), "\n")
  if( x$ws$homVar )
    cat("\nAlgorithm converged\n")
  else cat("\nAlgorithm converged in", x$ws$iter, 
           if(x$ws$iter>1) "iterations\n" else "iteration\n")
  invisible(x)
}

coef.nlreg <- function(object, ...)  
  object$coef

param <- function(object, ...)  
  UseMethod("param")

param.nlreg <- function(object, ...)  
  c(object$coef, object$varPar)

residuals.nlreg <- function(object, ...)  
  object$residuals

resid.nlreg <- residuals.nlreg

fitted.values.nlreg <- function(object, ...)
  object$fitted

fitted.nlreg <- fitted.values.nlreg

summary.nlreg <- function(object, observed = TRUE, 
                          correlation = FALSE, digits = NULL, ...)
{
  nlregObj <- object
  rc <- nlregObj$coef 
  vp <- nlregObj$varPar
  of <- nlregObj$offset
  if( nlregObj$ws$hoa )
  {
    md <- nlregObj$ws$md
    vd <- nlregObj$ws$vd
  }
  par <- nlregObj$ws$allPar
  mu <- nlregObj$fitted
  v <- nlregObj$weights
  .probl <- ( nlregObj$ws$homVar && !is.null(of) )
  if( .probl )
    .probl <- ( names(of) =="logs" )	
#  attach(nlregObj$data, warn.conflicts = FALSE)		## 20.05.13
  tmp <- as.list(par)
  tmp <- c(nlregObj$data, tmp)
  temp <- if( nlregObj$ws$hoa ) 
          {
          	formals(md) <- tmp
            do.call("md", tmp)
          }  
  	      else  
  	      {
  	        cat("\ndifferentiating mean function -- may take a while")
            tmp.fun <- Dmean(nlregObj)
          	formals(tmp.fun) <- tmp
            do.call("tmp.fun", tmp)
  	      }
  m1 <- attr(temp, "gradient") ; m1[!is.finite(m1)] <- 0
  m2 <- attr(temp, "hessian") ; m2[!is.finite(m2)] <- 0 
  if( !.probl )
  {
    tmp <- as.list(par)
    tmp <- c(nlregObj$data, tmp)
    temp <- if( nlregObj$ws$hoa )
            {
              formals(vd) <- tmp
              do.call("vd", tmp)
            }  
  	        else  
  	        {
  	          cat("\ndifferentiating variance function -- may take a while")
              tmp.fun <- Dvar(nlregObj)
              formals(tmp.fun) <- tmp
              do.call("tmp.fun", tmp)
            }
    v1 <- attr(temp, "gradient")  ; v1[!is.finite(v1)] <- 0    
    v2 <- attr(temp, "hessian") ; v2[!is.finite(v2)] <- 0 
  }
  else
  {
    v1 <- v2 <- NULL
  }
#  detach(nlreg.data)		## 20.05.13
  if( !nlregObj$ws$hoa )  cat("\n")
  info <- if( observed )  
            obsInfo.nlreg(nlregObj, par, mu, v, m1, m2, v1, v2)
 	    else  expInfo.nlreg(nlregObj, par, mu, v, m1, v1)
  s.e. <- sqrt(diag( (cov <- solve(qr(info))) ))
  dimnames(cov) <- list(names(c(rc,vp)), names(c(rc,vp)))
  idx1 <- 1:length(rc)
  if( !is.null(vp) )  
    idx2 <- length(rc) + 1:length(vp)	
  coef <- matrix(0, ncol = 4, nrow = length(rc))
  dimnames(coef) <- list(names(rc), 
                         c("Estimate", "Std. Error", "z value", 
                           "Pr(>|z|)"))
  coef[, 1] <- as.vector(rc)
  coef[, 2] <- s.e.[idx1]
  coef[, 3] <- coef[, 1]/coef[, 2]
  coef[, 4] <- 2 * pnorm( - abs(coef[, 3]) )
  if( !is.null(vp) )
  {
    varPar <- matrix(0, nrow = length(vp), ncol=2)
    dimnames(varPar) <- list(names(vp), c("Estimate", "Std. Error"))
    varPar[,1] <- vp
    varPar[,2] <- s.e.[idx2]	
  }
  else varPar <- NULL
  if( correlation )  
  {
    correl <- matvec( vecmat(1/s.e., cov),  1/s.e.)
    dimnames(correl) <- list(names(c(rc,vp)), names(c(rc,vp)))
  }
  else correl <- NULL
  summary <- list( coefficients = coef,
                   varPar       = varPar, 
                   offset       = nlregObj$offset,
                   residuals    = nlregObj$resid,
                   covariance   = cov,
                   correlation  = correl,
                   logLik       = nlregObj$logLik,
                   call         = nlregObj$call,
                   digits       = digits,
    		   ws           = nlregObj$ws )
  attr(summary, "class") <- c("summary.nlreg", "summary.nls")
  summary
}

print.summary.nlreg <- function(x, 
                         digits = max(3, getOption("digits")-3),
                         signif.stars = 
                                     getOption("show.signif.stars"),
                         quote = TRUE, ...)
{
  coef <- x$coef
  varPar <- x$varPar
  of <- x$offset
  correl <- x$correl
  if(missing(digits) && !is.null(x$digits))
     digits <- x$digits
  cat("Call: \n")
  dput(x$call)
  cat("\nRegression coefficients:\n")
  printCoefmat(coef, digits=digits, signif.stars=signif.stars)
  if( !is.null(varPar) )
  {
    cat("\nVariance parameters:\n")
    print(varPar, digits=digits)		
  }
  else cat("\nNo variance parameter\n")
  if( !is.null(of) )
  {
    cat("\nParameter of interest:\n")
    print(of, digits=digits)
  }
  else cat("\nNo interest parameter\n")	
  nobs <- length(x$residuals)      
  cat("\nTotal number of observations:", nobs)
  npar <- length(coef[,1]) + if( is.null(varPar) )  
                               0 else length(varPar[,1])
  cat("\nTotal number of parameters:", npar)
  cat("\n-2*Log Likelihood", format(-2 * x$logLik, digits=digits), 
      "\n")
  if( x$ws$homVar )
    cat("\nAlgorithm converged\n")
  else cat("\nAlgorithm converged in", x$ws$iter, 
           if(x$ws$iter>1) "iterations\n" else "iteration\n")
  if( !is.null(correl) ) 
  {
    p <- dim(correl)[2]
    if(p > 1) 
    {
      cat("\nCorrelation of parameters:\n")
      ll <- lower.tri(correl)
      correl[ll] <- format(round(correl[ll], digits))
      correl[!ll] <- ""
      print(correl[-1, - p, drop=FALSE], quote=FALSE, digits=digits)
    }
  }
  invisible(x)
}

var2cor <- function(object, ...)  
  UseMethod("var2cor")

var2cor.nlreg <- function(object, ...)
{
  object <- summary(object, corr=TRUE)
  var2cor(object)
}

var2cor.summary.nlreg <- function(object, ...)
{
  if( !is.null(object$correlation) )
    object$correlation
  else
  {
    cov <- object$covariance    
    s.e. <- sqrt(diag(cov))
    correl <- matvec( vecmat(1/s.e., cov),  1/s.e.)
    dimnames(correl) <- dimnames(object$covariance)
    correl
 }
}

nlreg.diag <- function(fitted, hoa = TRUE, infl = TRUE, trace = FALSE)
{
  nlregObj <- fitted
  rc <- nlregObj$coef 
  vp <- nlregObj$varPar
  of <- nlregObj$offset
  dupl <- nlregObj$data$dupl
  par <- nlregObj$ws$allPar
  mu <- nlregObj$fitted
  v <- nlregObj$weights
  .probl <- ( nlregObj$ws$homVar && !is.null(of) )
  if(.probl)
    .probl <- ( names(of) =="logs" )
  followUp <- trace
  wantHoa <- hoa
  if(wantHoa)
  {
    nlregHoa <- nlregObj
    nlregHoa$meanFun[[2]] <- call("+", nlregHoa$meanFun[[2]], 
                                  call("*", as.name("phi"), 
                                  as.name(".addVar")))
    nlregHoa$coef <- c(nlregHoa$coef, phi=0)
    nlregHoa$ws$allPar <- c(nlregHoa$ws$allPar, phi=0)
    cat("\ndifferentiating mean function -- may take a while")
    md <- Dmean(nlregHoa)
    if(!.probl)
    {
      cat("\ndifferentiating variance function -- may take a while")  
      vd <- Dvar(nlregHoa)
    }
    cut.idx <- length(rc) + 1		 
    .addVar <- rep(0, length(nlregObj$data$repl))  
    nlregHoa$data <- c(nlregHoa$data, list(.addVar=.addVar))
#    attach(nlregHoa$data, warn.conflicts = FALSE)		## 20.05.13
    new.par <- nlregHoa$ws$allPar
    tmp <- as.list(new.par)
    tmp <- c(nlregHoa$data, tmp)
    formals(md) <- tmp
    temp <- do.call("md", tmp)
    m1.0 <- attr(temp, "gradient") ; m1.0[!is.finite(m1.0)] <- 0
    m2.0 <- attr(temp, "hessian") ; m2.0[!is.finite(m2.0)] <- 0	
    m1.cut <- m1.0[,-cut.idx, drop=FALSE]
    m2.cut <- m2.0[,-cut.idx, -cut.idx, drop=FALSE]
    if(!.probl) 
    {
      formals(vd) <- tmp
      temp <- do.call("vd", tmp) 
      v1.0 <- attr(temp, "gradient") ; v1.0[!is.finite(v1.0)] <- 0
      v2.0 <- attr(temp, "hessian") ; v2.0[!is.finite(v2.0)] <- 0
      if( nlregObj$ws$xVar )
      {
	v1.cut <- v1.0[,-cut.idx, drop=FALSE]
	v2.cut <- v2.0[,-cut.idx, -cut.idx, drop=FALSE]
      }
      else  { v1.cut <- v1.0 ; v2.cut <- v2.0 } 
    }
    else 
    { 
      v1.0 <- v2.0 <- v1.cut <- v2.cut <- NULL 
    }
    M1.0 <- m1.cut[dupl,]
#    detach(nlregHoa$data)		## 20.05.13
    i.0 <- expInfo.nlreg(nlregObj, par, mu, v, m1.cut, v1.cut)
    j.0 <- obsInfo.nlreg(nlregObj, m1=m1.cut, m2=m2.cut, 
	                 v1=v1.cut, v2=v2.cut)
#    if( infl && !nlregObj$ws$homVar && !is.null(nlregObj$varPar ) )
    if( infl && !is.null(nlregObj$varPar ) )		## 14.06.13
    {
      mle.rc <- solve(j.0[names(rc), names(rc), drop=FALSE], 
                      j.0[names(rc), names(vp), drop=FALSE])
      mle.vp <- solve(j.0[names(vp), names(vp), drop=FALSE],   
                      j.0[names(vp), names(rc), drop=FALSE])
    }
    j.0 <- det(j.0)
    newCall <- match.call(nlreg, nlregObj$call)
    newCall$formula[[3]] <- nlregHoa$meanFun[[2]]
    newStart <- as.list(newCall$start)
    newStart <- c(newStart, list(phi=0))
    newStart <- as.call(newStart)
    newCall$start <- newStart
    newCall$trace <- FALSE
    if( nlregObj$ws$hoa )
      newCall$hoa <- FALSE
    newCall$data <- NULL
###    newData <- eval(nlregObj$ws$frame)		## 24.05.13
    if(!is.null(nlregObj$ws$frame))
      newData <- eval(nlregObj$ws$frame)
    else
    {
      .allVars <- all.vars(nlregObj$call$formula)
      .start <- names(nlregObj$call$start)
      .modVars <- .allVars[!match(.allVars, .start, nomatch=0)]	
      newData <- sapply(.modVars, function(x) eval(as.name(x)))
    }	    
    logLik.0 <- nlregObj$logLik
    phi.h <- logLik.1 <- j.1 <- i.1 <- S.1 <- q.1 <- rep(NA, length(v))
  }
  else
  {
#    if( infl && !nlregObj$ws$homVar )
    if( infl )		## 14.06.13
    {
      md <- if(nlregObj$ws$hoa )  
              nlregObj$ws$md
            else  
            {
              cat("\ndifferentiating mean function -- may take a while")
              Dmean(nlregObj)
            }
      vd <- if(nlregObj$ws$hoa )  
              nlregObj$ws$vd
            else  
            {
              cat("\ndifferentiating variance function -- may take a while")
              Dvar(nlregObj)
            }
    }
    else
    {
      md <- if(nlregObj$ws$hoa )  
              nlregObj$ws$md
            else  
            {
              cat("\ndifferentiating mean function -- may take a while")
              Dmean(nlregObj, hessian=FALSE)
            }
      if( !.probl )
	  vd <- if(nlregObj$ws$hoa )  
                nlregObj$ws$vd
              else  
              {
                cat("\ndifferentiating variance function -- may take a while")
                Dvar(nlregObj, hessian=FALSE)
              }
    }
#    attach(nlregObj$data, warn.conflicts = FALSE)		## 20.05.13
    tmp <- as.list(nlregObj$ws$allPar)
    tmp <- c(nlregObj$data, tmp)
    formals(md) <- tmp
    temp <- do.call("md", tmp)
    m1.0 <- attr(temp, "gradient") ; m1.0[!is.finite(m1.0)] <- 0
#    if( infl && !nlregObj$ws$homVar )
    if( infl )		## 14.06.13
    { 
      m2.0 <- attr(temp, "hessian") ; m2.0[!is.finite(m2.0)] <- 0 
    }
    if(!.probl) 
    {
      formals(vd) <- tmp
      temp <- do.call("vd", tmp)
      v1.0 <- attr(temp, "gradient") ; v1.0[!is.finite(v1.0)] <- 0
#      if( infl && !nlregObj$ws$homVar )
      if( infl )		## 14.06.13
      { 
        v2.0 <- attr(temp, "hessian") ; v2.0[!is.finite(v2.0)] <- 0 
      }
    }
    else 
    { 
      v1.0 <- v2.0 <- NULL 
    }
    M1.0 <- m1.0[dupl,]
#    detach(nlregObj$data)		## 20.05.13
    i.0 <- expInfo.nlreg(nlregObj, par, mu, v, m1.0, v1.0)
#    if( infl && !nlregObj$ws$homVar && !is.null(nlregObj$varPar) )
    if( infl && !is.null(nlregObj$varPar) )		## 14.06.13
    { 
      j.0 <- obsInfo.nlreg(nlregObj, m1=m1.0, m2=m2.0, v1=v1.0, 
                           v2=v2.0)
      mle.rc <- solve(j.0[names(rc), names(rc), drop=FALSE], 
                      j.0[names(rc), names(vp), drop=FALSE])
      mle.vp <- solve(j.0[names(vp), names(vp), drop=FALSE], 
                      j.0[names(vp), names(rc), drop=FALSE])
    }
  }
#  if( infl && !nlregObj$ws$homVar )
  if( infl )		## 14.06.13
  {
    newFrame <- eval(nlregObj$ws$frame)
    newFrame <- c(newFrame, as.list(nlregObj$coef))
    newFrame <- c(newFrame, as.list(nlregObj$varPar))
    newFrame <- c(newFrame, as.list(nlregObj$offset))
    logLik1 <- nlregObj$logLik
    ld <- ld.rc <- ld.vp <- rep(0, length(v))
  }
  res <- nlregObj$resid	                    ## standardized residuals
  s.i.0 <- solve(qr(i.0))
  dimnames(s.i.0) <- dimnames(i.0)
  hh <- diag( M1.0 %*% s.i.0[names(rc), names(rc), drop=FALSE] %*% 
                       t(M1.0) ) / v	                 ## leverages
  ha <- diag( M1.0 %*% solve( crossprod( M1.0, vecmat(1/v, M1.0) ),
                              t(M1.0) ) ) / v 
                                              ## approximate leverages
  r.p <- nlregObj$residuals / sqrt( 1 - hh )	
                                     ## generalized Pearson residuals
  r.s <- nlregObj$residuals / sqrt( 1 - ha )	
                                 ## approximate studentized residuals
  cook <- hh/(length(rc)*(1-hh))*r.p^2		
                 ## approximate influence for regression coefficients
  if(followUp)  
    cat("\n")
  for( sset in 1:length(v) )
  {
    if(wantHoa)
    {
      .addVar <- diag(1, length(v))[, sset]
      ..newData2 <- data.frame(cbind(newData, .addVar))
      lv <- lm((nlregObj$residual*sqrt(nlregObj$weights))~.addVar-1)$coef
      newCall$start["phi"] <- signif(lv, 4)
###      attach(..newData2, warn.conflicts = FALSE)		## 24.05.13
      newCall$data <- as.name("..newData2")
       ctrl <- try(eval(newCall), silent=TRUE)
      if( class(ctrl)[1] != "try-error" )
      {
        nlreg.temp <- eval(newCall)
        lv <- coef(nlreg.temp)["phi"]
        phi.h[sset] <- nlreg.temp$coef["phi"]
        logLik.1[sset] <- nlreg.temp$logLik
#        attach(nlreg.temp$data, warn.conflicts = FALSE)		## 20.05.13
        tmp <- as.list(nlreg.temp$ws$allPar)
        tmp <- c(nlreg.temp$data, list(.addVar=.addVar), tmp)
###        tmp <- c(nlreg.temp$data, tmp)
        formals(md) <- tmp
        temp <- do.call("md", tmp)
        m1.1 <- attr(temp, "gradient") ; m1.1[!is.finite(m1.1)] <- 0
        m2.1 <- attr(temp, "hessian") ; m2.1[!is.finite(m2.1)] <- 0
        tmp <- as.list(nlregHoa$ws$allPar)
        tmp <- c(nlreg.temp$data, list(.addVar=.addVar), tmp)
###        tmp <- c(nlreg.temp$data, tmp)
        formals(md) <- tmp
        temp <- do.call("md", tmp)
        m1.0 <- attr(temp, "gradient") ; m1.0[!is.finite(m1.0)] <- 0
        if(!.probl) 
        {
          tmp <- as.list(nlreg.temp$ws$allPar)
          tmp <- c(nlreg.temp$data, list(.addVar=.addVar), tmp)
###          tmp <- c(nlreg.temp$data, tmp)
          formals(vd) <- tmp
          temp <- do.call("vd", tmp)
          v1.1 <- attr(temp, "gradient") ; v1.1[!is.finite(v1.1)] <- 0
          v2.1 <- attr(temp, "hessian") ; v2.1[!is.finite(v2.1)] <- 0
          tmp <- as.list(nlregHoa$ws$allPar)
          tmp <- c(nlreg.temp$data, list(.addVar=.addVar), tmp)
###          tmp <- c(nlreg.temp$data, tmp)
          formals(vd) <- tmp
          temp <- do.call("vd", tmp)
	      v1.0 <- attr(temp, "gradient") ; v1.0[!is.finite(v1.0)] <- 0
        }
        else v1.1 <- v2.1 <- v1.0 <- NULL
#        detach(nlreg.temp$data)		## 20.05.13
###        detach(..newData2)		## 24.05.13
        j.1[sset] <- det(obsInfo.nlreg(nlreg.temp, m1=m1.1, m2=m2.1, 
  	                               v1=v1.1, v2=v2.1))
        i.1[sset] <- det(expInfo.nlreg(nlreg.temp, m1=m1.1, v1=v1.1))
        S.temp <- Shat.nlreg(nlreg.temp, nlregObj, par.0=new.par,  
	                     m1.1=m1.1, m1.0=m1.0, v1.1=v1.1, v1.0=v1.0)
        q.temp <- qhat.nlreg(nlreg.temp, nlregObj, par.0=new.par, 
                             m1.1=m1.1, v1.1=v1.1)
        s.S.temp <- solve(qr(S.temp))
        dimnames(s.S.temp) <- dimnames(S.temp)
        q.1[sset] <- s.S.temp["phi",] %*% q.temp
        S.1[sset] <- det(S.temp)
      }
    }
#    if( infl && !nlregObj$ws$homVar )		## 14.06.13
    if( infl )
    {
      nlreg.temp <- update(nlregObj, subset=c(-sset), trace=FALSE)
      newFrame[names(nlreg.temp$coef)] <- nlreg.temp$coef
      newFrame[names(nlreg.temp$varPar)] <- nlreg.temp$varPar
      num <- eval( call("-", nlregObj$call$formula[[2]], 
                        nlregObj$call$formula[[3]]), envir=newFrame )
      num <- num^2
      den <- eval( nlregObj$varFun[[2]], envir=newFrame )
      logLik0 <- -( sum(log(den)) + sum(num/den) + 
                    length(num)*log(2*pi) )/2
      ld[sset] <- 2/length(c(rc,vp))*(logLik1-logLik0)
      if( !is.null(nlregObj$varPar) )
      {
        rc.cond <- rc + mle.rc %*% (vp - nlreg.temp$varPar)
        vp.cond <- vp + mle.vp %*% (rc - nlreg.temp$coef)	
        newFrame[names(nlreg.temp$varPar)] <- vp.cond
        num <- eval( call("-", nlregObj$call$formula[[2]], 
                          nlregObj$call$formula[[3]]), 
                     envir=newFrame )
        num <- num^2
        den <- eval( nlregObj$varFun[[2]], envir=newFrame )
        logLik0 <- -( sum(log(den)) + sum(num/den) + 
                      length(num)*log(2*pi) )/2
        ld.rc[sset] <- 2/length(rc)*(logLik1-logLik0)
        newFrame[names(nlreg.temp$coef)] <- rc.cond
        newFrame[names(nlreg.temp$varPar)] <- nlreg.temp$varPar
        num <- eval( call("-", nlregObj$call$formula[[2]], 
                          nlregObj$call$formula[[3]]), 
                     envir=newFrame )
        num <- num^2
        den <- eval( nlregObj$varFun[[2]], envir=newFrame )
        logLik0 <- -( sum(log(den)) + sum(num/den) +
                      length(num)*log(2*pi) )/2
        ld.vp[sset] <- 2/length(vp)*(logLik1-logLik0)
      }
    }
    if(followUp)
    {
#      if(wantHoa || (infl && !nlregObj$ws$homVar))
      if(wantHoa || infl)		## 14.06.13
	cat(paste("\niteration", sset, "out of", length(v)))
      if( wantHoa )
      {
        if( class(ctrl)[1] == "try-error" )
          cat("\n\'nlreg\' did not converge: NAs produced in diagnotics")
      }
    }
    else
      if( wantHoa )
      {
        if( class(ctrl)[1] == "try-error" )
          cat(paste("\n\'nlreg\' did not converge in iteration",
                    sset, "out of", length(v), 
                    ": NAs produced in diagnotics \n"))
      }
  }
  if(followUp)  
#    if(wantHoa || (infl && !nlregObj$ws$homVar))  cat("\n")    
	 if(wantHoa || infl)  cat("\n")		## 14.06.13
  if(wantHoa)  
  {
    r.j <- sign(phi.h) * sqrt(2*(logLik.1-logLik.0))	
                                                 ## deletion residuals
    rs.j <- r.j - log(abs(i.1*sqrt(j.0)*r.j/
                                      sqrt(j.1)/S.1/q.1))/r.j
                                                 ## r*-type residuals
  }
  diagObj <- list(fitted = nlregObj$fitted, 
                   resid = res, rp = r.p, rs = r.s)
  if(wantHoa)
    diagObj <- c(diagObj, list( rj = r.j, rsj = rs.j ) )
  diagObj <- c(diagObj, list( h = hh, ha = ha, cook = cook ) )
#  if( infl && !nlregObj$ws$homVar )
  if( infl )		## 14.06.13
  {
    diagObj <- c(diagObj, list( ld = ld) )
    if( !is.null(nlregObj$varPar) )
      diagObj <- c(diagObj, list( ld.rc = ld.rc, ld.vp = ld.vp ))
  }
  diagObj <- c( diagObj, npar=length(rc) )
  class(diagObj) <- "nlreg.diag"
  diagObj
}	

nlreg.diag.plots <- function(fitted, which = "all", subset = NULL, 
                             iden = FALSE, labels = NULL, hoa = TRUE, 
                             infl = TRUE, trace = FALSE, ret = FALSE, 
                             ...)
{
  if( any(match("nlreg.diag", class(fitted), nomatch=FALSE)) )
  {
    plot(fitted, which=which, subset=subset, iden=iden, labels=labels,
         ...)
    if(ret) fitted
  }
  else	if( any(match("nlreg", class(fitted), nomatch=FALSE)) )
	{
	  objNlreg <- nlreg.diag(fitted, hoa=hoa, infl=infl, 
                                  trace=trace) 
	  plot(objNlreg, which=which, subset=subset, iden=iden, 
               labels=labels, ...)
	  if (ret) objNlreg
	}
	else
  	  stop("missing \"nlreg\" or \"nlreg.diag\" object")
}	

plot.nlreg.diag <- function(x, which = "all", subset = NULL, 
                            iden = FALSE, labels = NULL, ...)
{
  nlregdiag <- x
  if( attr(nlregdiag, "class") != "nlreg.diag" )
    stop("\nmissing \"nlreg.diag\" object")
  if( is.null(subset) )
    subset <- c( 1:length(nlregdiag$h) )
  else if( is.logical(subset) )
	 subset <- ( 1:length(subset) )[subset]
       else if( is.numeric(subset) && all(subset < 0) )
	      subset <- ( 1:(length(subset) + 
                             length(nlregdiag$h)) )[subset]
	    else if(is.character(subset)) 
	         {
		   if( is.null(labels) )
		     labels <- subset
		   subset <- seq(along = subset)
	         }
  choices <- c("Summary", 
	       "Studentized residuals against fitted values",
	       "r* residuals against fitted values", 
               "Normal QQ-plot of studentized residuals",
	       "Normal QQ-plot of r* residuals", 
               "Cook statistic against h/(1-h)",
               "Global influence against h/(1-h)",		
	       "Cook statistic against observation number",
	       "Influence measures against observation number\n")
  tmenu <- paste("plot:", choices)
  if( which == "all" )
    pick <- menu(tmenu, 
                 title = "\n Make a plot selection (or 0 to exit)\n")
  else if( !match(which, 2:length(choices), nomatch=FALSE) )
	 stop("choice not valid") 
       else  pick <- which
  if( pick == 0 )
    stop("no graph required ! ")
  old.par <- par(no.readonly = TRUE) ; on.exit(par(old.par))
  repeat
  {
    switch(pick, 
           "1" = { par(pty="s", mfrow=c(1,1))
  		   close.screen(all.screens = TRUE)
		   split.screen(figs = c(2, 2))
	           screen(1)       
##  Plot the studentized residuals against the fitted values
                   x1 <- nlregdiag$fitted
                   plot(x1, nlregdiag$rp, xlab = "Fitted values", 
                        ylab = "Studentized residuals", ...)
		   screen(2)       
## Plot a normal QQ-plot of either the studentized residuals or the r*
## residuals
                   y2 <- if( !is.null(nlregdiag$rsj) )  nlregdiag$rsj
			 else  nlregdiag$rp
                   x2 <- qnorm(ppoints(length(y2)))[rank(y2)]
                   plot(x2, y2, 
                        xlab = paste("Quantiles of standard normal"), 
                        ylab = paste("Ordered", 
                                     if( !is.null(nlregdiag$rsj) ) 
                                       "r* residuals" 
                                     else "studentized residuals"), 
                        xlim=c(-3,3), ylim=c(-3,3), ...)
                   abline(0, 1, lty = 2)
		   screen(3)
## Plot either the Cook statistics or the global influence measure
## against h/(1-h) and draw line to highlight possible influential and
## high leverage points.
                   hh <- nlregdiag$h/(1 - nlregdiag$h)
		   y3 <- if( is.null(nlregdiag$ld) )  
                           nlregdiag$cook else nlregdiag$ld
                   plot(hh, y3, xlab = "h/(1-h)",  
	                ylab = if( is.null(nlregdiag$ld) ) 
                                 "Cook statistic" 
                               else "Global influence", ...)
	           rx <- range(hh)
                   ry <- range(y3)
                   npar <- nlregdiag$npar
                   nobs <- length(nlregdiag$h)
                   cooky <- 8/(nobs - 2 * npar)
                   hy <- (2 * npar)/(nobs - 2 * npar)
                   if((cooky >= ry[1]) && (cooky <= ry[2]))
                     abline(h = cooky, lty = 2)
                   if((hy >= rx[1]) && (hy <= rx[2]))
                     abline(v = hy, lty = 2)
                   screen(4)
## Plot either the Cook statistics or the influence measures against
## the observation number in the original data set.
		   y4 <- if( is.null(nlregdiag$ld) )  nlregdiag$cook 
                         else nlregdiag$ld
		   plot(subset, y4, xlab = "Case", type="h",
			ylab = if( is.null(nlregdiag$ld) ) 
                                 "Cook statistic" 
                               else "Global influence", ...)
		   if((cooky >= ry[1]) && (cooky <= ry[2]))
		     abline(h = cooky, lty = 2)
		   xx <- list(x1, x2, hh, subset)
		   yy <- list(nlregdiag$rp, y2, y3, y4)
                   if(is.null(labels))
                     labels <- if( !is.null(names(nlregdiag$fitted)) )
				 names(nlregdiag$fitted) 
			       else paste(1:length(nlregdiag$h))
 	           yes <- iden
		   while(yes) 
		   {
##  If interaction with the plots is required then ask the user which
##  plot they wish to interact with and then run identify() on that
##  plot.  When the user terminates identify(), reprompt until no
##  further interaction is required and the user inputs a 0.
                     cat("****************************************************\n")
 	             cat("Please Input a screen number (1,2,3 or 4)\n")
                     cat("0 will terminate the function \n")
                     num <- scan(n = 1)
		     check <- ((length(num) > 0) && ((num == 1) || (
                               num == 2) || (num == 3) || (num == 4)))
		     if(check) 
		     {
                       cat(paste("Interactive Identification for screen", 
                                 num, "\n"))
                       cat("left button = Identify, center button = Exit\n")
		       screen(num, new = FALSE)
		       identify(xx[[num]], yy[[num]], labels, ...)
		     }
		     else yes <- FALSE
		   }
		   close.screen(all.screens=TRUE)
		   par(ask=FALSE)
		 },

           "2" = { par(pty="s", mfrow=c(1,1))
##  Plot the studentized residuals against the fitted values
	           x1 <- nlregdiag$fitted
                   plot(x1, nlregdiag$rp, xlab = "Fitted values", 
	                ylab = "Studentized residuals", ...)
                   xx <- list(x1)
                   yy <- list(nlregdiag$rp)
		 },
	   "3" = { par(pty="s", mfrow=c(1,1))
##  Plot the r* residuals (if available) against the fitted values
		   if( is.null(nlregdiag$rsj) )
		   {
		     plot(0,0, type="n", axes=FALSE, xlab="", ylab="")
		     text(0,0, "Not available")         
		   }
		   else
		   {
                     x2 <- nlregdiag$fitted
                     plot(x2, nlregdiag$rsj, xlab = "Fitted values", 
                          ylab = "r* residuals", ...)
                     xx <- list(x2)
                     yy <- list(nlregdiag$rsj)
		   }
	         },
           "4" = { par(pty="s", mfrow=c(1,1))
##  Plot a normal QQ-plot of the studentized residuals
                   y3 <- nlregdiag$rp
                   x3 <- qnorm(ppoints(length(y3)))[rank(y3)]
                   plot(x3, y3, xlab = "Quantiles of standard normal",
                        ylab = "Ordered studentized residuals", 
                        xlim=c(-3,3), ylim=c(-3,3), ...)
                   abline(0, 1, lty = 2)
                   xx <- list(x3)
                   yy <- list(y3)
                 },
           "5" = { par(pty="s", mfrow=c(1,1))
##  Plot a normal QQ-plot of the r* residuals (if available)
		   if( is.null(nlregdiag$rsj) )
		   {
		     plot(0,0, type="n", axes=FALSE, xlab="", ylab="")
                     text(0,0, "Not available")
		   }
		   else
	           {                        
                     y4 <- nlregdiag$rsj
                     x4 <- qnorm(ppoints(length(y4)))[rank(y4)]
                     plot(x4, y4, 
                          xlab = "Quantiles of standard normal", 
                          ylab = "Ordered r* residuals", ylim=c(-3,3),
                          xlim=c(-3,3), ...) 
                     abline(0, 1, lty = 2)
                     xx <- list(x4)
                     yy <- list(y4)
	           }
                 },
           "6" = { par(pty="s", mfrow=c(1,1))
## Plot the Cook statistics against h/(1-h) and draw line to highlight
## possible influential and high leverage points.
                   hh <- nlregdiag$h/(1 - nlregdiag$h)
		   y5 <- nlregdiag$cook 
	           plot(hh, y5, xlab = "h/(1-h)", 
	                ylab = "Cook statistic", ...)
	           rx <- range(hh)
                   ry <- range(y5)
                   npar <- nlregdiag$npar
                   nobs <- length(nlregdiag$h)
                   cooky <- 8/(nobs - 2 * npar)
                   hy <- (2 * npar)/(nobs - 2 * npar)
                   if((cooky >= ry[1]) && (cooky <= ry[2]))
                     abline(h = cooky, lty = 2)
                   if((hy >= rx[1]) && (hy <= rx[2]))
                     abline(v = hy, lty = 2)
                   xx <- list(hh)
                   yy <- list(nlregdiag$cook)
                 },
           "7" = { par(pty="s", mfrow=c(1,1))
## Plot the global influence (if available) against h/(1-h) and draw
## line to highlight possible influential and high leverage points.
		   if( is.null(nlregdiag$ld) )
		   {
		     plot(0,0, type="n", axes=FALSE, xlab="", ylab="")
		     text(0,0, "Not available")
		   }	
		   else
		   {
		     hh <- nlregdiag$h/(1 - nlregdiag$h)
 		     y6 <- nlregdiag$ld
	             plot(hh, y6, xlab = "h/(1-h)", 
	                  ylab = "Global influence", ...)
  	             rx <- range(hh)
                     ry <- range(y6)
                     npar <- nlregdiag$npar
                     nobs <- length(nlregdiag$h)
                     cooky <- 8/(nobs - 2 * npar)
                     hy <- (2 * npar)/(nobs - 2 * npar)
                     if((cooky >= ry[1]) && (cooky <= ry[2]))
                       abline(h = cooky, lty = 2)
                     if((hy >= rx[1]) && (hy <= rx[2]))
                       abline(v = hy, lty = 2)
                     xx <- list(hh)
                     yy <- list(nlregdiag$ld)
	           }
                 },
	   "8" = { par(pty="s", mfrow=c(1,1))
## Plot the Cook statistics against the observation number in the
## original data set.
		   plot(subset, nlregdiag$cook, xlab = "Case", 
                        ylab = "Cook statistic", type="h", ...)
                   ry <- range(nlregdiag$cook)
                   npar <- nlregdiag$npar
                   nobs <- length(nlregdiag$h)
                   cooky <- 8/(nobs - 2 * npar)
 	           if((cooky >= ry[1]) && (cooky <= ry[2]))
                   xx <- list(subset)
                   yy <- list(nlregdiag$cook)
                 },
	   "9" = { par(pty="s", mfrow=c(1,1))
## Plot the influence measures (if available) against the observation
## number in the original data set.
		   if( is.null(nlregdiag$ld) )
		   {
		     plot(0,0, type="n", axes=FALSE, xlab="", ylab="")
		     text(0,0, "Not available")
		   }	
		   else
		   {
	             if( is.null(nlregdiag$ld.rc) )
                     {
		       plot(subset, nlregdiag$ld, xlab = "Case", 
                            ylab = "", main = "Global Influence", 
                            type="h", ...)
                       ry <- range(nlregdiag$ld)
		       npar <- nlregdiag$npar
		       nobs <- length(nlregdiag$h)
	               cooky <- 8/(nobs - 2 * npar)
		       if((cooky >= ry[1]) && (cooky <= ry[2]))
			 abline(h = cooky, lty = 2)
		       if(is.null(labels))
			 labels <- if( !is.null(names(nlregdiag$fitted)) )
			              names(nlregdiag$fitted) 
				   else paste(1:length(nlregdiag$h))
##  If interaction with the plots is required then ask the user which
##  plot they wish to interact with and then run identify() on that
##  plot.  When the user terminates identify(), reprompt until no
##  further interaction is required and the user inputs a 0.
		       cat("****************************************************\n")
		       cat("Interactive Identification\n")
	       cat("left button = Identify, center button = Exit\n")
		       identify(subset, nlregdiag$ld, labels, ...)
		     }
		     else
	             {			
		       split.screen(figs=c(1,3))
		       screen(1)
		       plot(subset, nlregdiag$ld, xlab = "Case", 
                            ylab = "", main = "Global Influence", 
                            type="h", ...)
                       ry <- range(nlregdiag$ld)
	               npar <- nlregdiag$npar
		       nobs <- length(nlregdiag$h)
	               cooky <- 8/(nobs - 2 * npar)
		       if((cooky >= ry[1]) && (cooky <= ry[2]))
			 abline(h = cooky, lty = 2)
 		       screen(2)
		       plot(subset, nlregdiag$ld.rc, xlab = "Case", 
                            ylab = "",
                         main = "Partial Influence\n(regr. coeff.)\n",
                            type="h", ...)
		       screen(3)
		       plot(subset, nlregdiag$ld.vp, xlab = "Case", 
                            ylab = "",
                          main = "Partial Influence\n(var. param.)\n",
                            type="h", ...)
		       xx <- list(subset, subset, subset)
		       yy <- list(nlregdiag$ld, nlregdiag$ld.rc, 
                                  nlregdiag$ld.vp)
		       if(is.null(labels))
		 labels <- if( !is.null(names(nlregdiag$fitted)) )
				     names(nlregdiag$fitted) 
				   else paste(1:length(nlregdiag$h))
		       yes <- iden
		       while(yes) 
		       {
##  If interaction with the plots is required then ask the user which
##  plot they wish to interact with and then run identify() on that
##  plot.  When the user terminates identify(), reprompt until no
##  further interaction is required and the user inputs a 0.
	                 cat("****************************************************\n")
		 cat("Please Input a screen number (1,2 or 3)\n")
			 cat("0 will terminate the function \n")
			 num <- scan(n = 1)
		 check <- ((length(num) > 0) && ((num == 1) || 
                                   (num == 2) || (num == 3)))
 		         if(check) 
		         {
		   cat(paste("Interactive Identification for screen", 
                             num, "\n"))
	   cat("left button = Identify, center button = Exit\n")
 			   screen(num, new = FALSE)
			   identify(xx[[num]], yy[[num]], labels, ...)
		         }
		         else yes <- FALSE
		       }	
		       close.screen(all.screens=TRUE)
		       par(ask=FALSE)
                     }
		   }} )     
    if( (pick != 1) && (pick != 9) && iden )
    {           
      if( !( (((pick==3) || (pick==5)) && is.null(nlregdiag$rsj)) || 
              ((pick==7) && is.null(nlregdiag$ld)) ) )
      { 
	if(is.null(labels))
	  labels <- if( !is.null(names(nlregdiag$fitted)) )
		      names(nlregdiag$fitted) 
		    else paste(1:length(nlregdiag$h))
##  If interaction with the plots is required then ask the user which
##  plot they wish to interact with and then run identify() on that
##  plot.  When the user terminates identify(), reprompt until no
##  further interaction is required and the user inputs a 0.
	cat("****************************************************\n")
	cat("Interactive Identification\n")
	cat("left button = Identify, center button = Exit\n")
	identify(xx[[1]], yy[[1]], labels, ...)
      }
    }
    if( which == "all" )
      pick <- menu(tmenu, 
                  title = "\n Make a plot selection (or 0 to exit)\n")
    if( (pick == 0) || (which != "all") )
    {
      invisible(close.screen(all.screens=TRUE))    
      break
    }
  }       
  invisible()
}

profile.nlreg <- function(fitted, offset = "all", hoa = TRUE,       
                          precision = 6, signif = 30, n = 50, 
                          omit = 0.5, trace = FALSE, md, vd, 
                          all = FALSE, ...)
{
  m <- match.call()
  nlregObj <- fitted 
  if( !is.null(nlregObj$offset) )
    stop("\noffset parameter is already fixed!")
  offsetName <- paste(substitute(offset))
  if(offsetName == "all")
    all.profiles.nlreg(fitted=nlregObj, hoa=hoa, signif=signif, 
                       n=n, precision=precision, omit=omit, 
                       trace=trace, call=m, ...)
  else
  {
    if( !match(offsetName, names(nlregObj$ws$allPar), nomatch=FALSE) )
      stop("\noffset parameter does not appear in formula")
    signif <- 2*as.integer(signif/2)		
    followUp <- trace				
    wantHoa <- hoa		
    arg.list <- match.call(nlreg, nlregObj$call)
    arg.list$trace <- FALSE
    new.start <- nlregObj$ws$allPar[names(nlregObj$ws$allPar)!=
                                    "logs"] 
    if( offsetName != "logs" )
    {          
      idx <- match( offsetName, names(new.start) )
      new.start <- new.start[-idx]
    }
    if( nlregObj$ws$homVar )
      arg.list$weights <- NULL
    rc <- if( (idx1 <- match(offsetName, names(nlregObj$coef), 
                             nomatch=0)) )
 	    nlregObj$coef[-idx1]
          else nlregObj$coef	
    vp <- if( (idx2 <- match(offsetName, names(nlregObj$varPar), 
                             nomatch=0)) )
            nlregObj$varPar[-idx2]
          else nlregObj$varPar
    smaller <- c(names(rc), names(vp))
    if( missing(md) )
      if( nlregObj$ws$hoa )  
        md <- nlregObj$ws$md
      else  
      {
        cat("\ndifferentiating mean function -- may take a while")
        md <- Dmean(nlregObj)
      }
    if( missing(vd) )
      if( nlregObj$ws$hoa )  
        vd <- nlregObj$ws$vd
      else  
      {
    	cat("\ndifferentiating variance function -- may take a while")
        vd <- Dvar(nlregObj)  
      }
#    attach(nlregObj$data, warn.conflicts = FALSE)		## 20.05.13
    tmp <- as.list(nlregObj$ws$allPar)
    tmp <- c(nlregObj$data, tmp)
    formals(md) <- tmp
    temp <- do.call("md", tmp)
    m1.1 <- attr(temp, "gradient") ; m1.1[!is.finite(m1.1)] <- 0
    m2.1 <- attr(temp, "hessian") ; m2.1[!is.finite(m2.1)] <- 0
    formals(vd) <- tmp
    temp <- do.call("vd", tmp)
    v1.1 <- attr(temp, "gradient") ; v1.1[!is.finite(v1.1)] <- 0 
    v2.1 <- attr(temp, "hessian") ; v2.1[!is.finite(v2.1)] <- 0
#    detach(nlregObj$data)		## 20.05.13
    logLik1 <- nlregObj$logLik
    mle <- if( idx1 )  nlregObj$coef[idx1]
             else nlregObj$varPar[idx2]		
    jh <- obsInfo.nlreg(nlregObj, m1=m1.1, m2=m2.1, v1=v1.1, v2=v2.1)
    if(all)  jh.all <- jh
    s.e. <- solve(qr(jh))
    dimnames(s.e.) <- dimnames(jh)
    s.e. <- sqrt(s.e.[offsetName, offsetName])
    if(wantHoa)
    {
      jh <- sqrt(det(jh))
      jh.t <- jh * s.e.
      ih <- expInfo.nlreg(nlregObj, m1=m1.1, v1=v1.1)
      ih.t <- det(ih[smaller, smaller, drop=FALSE])
      ih <- det(ih)
      rh <- resid(nlregObj)  
      wh <- sqrt(nlregObj$weights)      
      db <- ncol(m1.1)
      dp <- ncol(v1.1) + ifelse(!nlregObj$ws$xVar, db, 0) 
      npts <- nrow(m1.1)
      Dvh <- matrix(0, nrow=npts, ncol=dp) 
      if(nlregObj$ws$xVar)
        Dvh <- v1.1
      else Dvh[,(db+1):dp] <- v1.1
      Dmh <- m1.1
      if( !nlregObj$ws$missingData )
      {
        Dvh <- Dvh[nlregObj$data$dupl,]
        Dmh <- Dmh[nlregObj$data$dupl,]
      }	
      Dvh <- Dvh/2/wh  
      Dmvh <- t(Dvh*rh) 
      Dmvh[1:db,] <- Dmvh[1:db,] + t(Dmh)
      dimnames(Dmvh) <- list(names(nlregObj$ws$allPar), rep("", ncol(Dmvh)))
      lvh <- - matrix( rowSums( t(Dvh / wh) ), ncol=1)
      lthvh <- Dmvh %*% (aa <- cbind( vecmat( 1/wh^2, Dmh ) + 
                                        vecmat( rh/wh^2, Dvh[,1:db] ),
                                      vecmat( 2*rh/wh^2, Dvh[,(db+1):dp] ) ))
#      lthvh <- Dmvh %*% (aa <- cbind( vecmat( 1/wh^2, (Dmh + Dvh[,1:db]) ),
#                                      vecmat( 2*rh/wh^2, Dvh[,(db+1):dp] ) ))
      dimnames(lthvh) <- list(names(nlregObj$ws$allPar), 
                              names(nlregObj$ws$allPar))
      dimnames(aa) <- list(rep("", nrow(aa)), names(nlregObj$ws$allPar))
      lchvh <-  Dmvh[smaller,,drop=FALSE] %*% aa[,smaller,drop=FALSE]
      if( nlregObj$ws$homVar )
        lthvh[dp,1:db] <- lthvh[1:db,dp] <- 0
    }
    offsetVals <- c(mle-precision*s.e., mle+precision*s.e.)
    offsetVals <- c(offsetVals[1], 
                    offsetVals[1] + cumsum(rep(diff(offsetVals)/
                                           (signif-1),(signif-1))))
    offsetVals <- signif(offsetVals, 6)
    x <- r <- w <- rep(NA, signif)
    if(wantHoa)
    {
      rs.sk <- q.sk <- l.sk <- inf.sk <- np.sk <- rep(NA, signif)
      rs.fr <- q.fr <- l.fr <- inf.fr <- np.fr <- rep(NA, signif)
    }
    mat <- matrix( rep(NA, signif*(length(nlregObj$ws$allPar) +
                           ifelse(wantHoa, 4, 2))), 
                   nrow=length(offsetVals) )
    new.start.0 <- new.start 
    go.on <- TRUE ; next.idx <- signif/2
    counter <- 0 ; first.time <- TRUE ; limit <- 0
    if(followUp)  
      if(nlregObj$ws$hoa)  cat("\n")
    else  cat("\n\n")

    for(i in 1:2)
    {
      while( go.on && (next.idx!=limit) )
      {
  	of <- offsetVals[next.idx]
  	x[next.idx] <- mat[next.idx,1] <- of
  	names(of) <- offsetName
        arg.list$offset <- 
          eval(parse(text=paste("call(\"c\",", offsetName, "=", of, ")")))
        arg.list$start <- as.name("new.start.0")
  	if( nlregObj$ws$hoa )
  	  arg.list$hoa <- FALSE
        ctrl <- try(eval(arg.list), silent=TRUE)
        if( class(ctrl)[1] != "try-error" )
        {
          nlregObj0 <- eval(arg.list)
          logLik0 <- nlregObj0$logLik
          w[next.idx] <- (mle-of)/s.e.
          r[next.idx] <- sqrt( 2* (logLik1-logLik0) ) * sign(mle-of)
          if(wantHoa)
          {
#            attach(nlregObj0$data, warn.conflicts = FALSE)		## 20.05.13
            tmp <- as.list(nlregObj0$ws$allPar)
            tmp <- c(nlregObj0$data, tmp)
            formals(md) <- tmp
            temp <- do.call("md", tmp)
  	        m1.0 <- attr(temp, "gradient") ; m1.0[!is.finite(m1.0)] <- 0
            m2.0 <- attr(temp, "hessian") ; m2.0[!is.finite(m2.0)] <- 0
            m1.00 <- m1.0[,names(rc), drop=FALSE]
  	        m2.00 <- m2.0[,names(rc),names(rc), drop=FALSE]
  	        formals(vd) <- tmp
            temp <- do.call("vd", tmp)
  	        v1.0 <- attr(temp, "gradient") ; v1.0[!is.finite(v1.0)] <- 0
  	        v2.0 <- attr(temp, "hessian") ; v2.0[!is.finite(v2.0)] <- 0
  	        if( (offsetName == "logs") && nlregObj0$ws$homVar )
            { v1.00 <- v1.0 ; v2.00 <- v2.0 }
  	        else if( !nlregObj0$ws$xVar)
  	        {
  		      v1.00 <- v1.0[, names(vp), drop=FALSE]
  		      v2.00 <- v2.0[, names(vp), names(vp), drop=FALSE]
  	        }
  	        else
  	        {
  	          v1.00 <- v1.0[, c(names(rc), names(vp)), drop=FALSE]
  		      v2.00 <- v2.0[, c(names(rc), names(vp)), 
  	                          c(names(rc), names(vp)), drop=FALSE]
  	        }
#  	    detach(nlregObj0$data)		## 20.05.13
            jt <- obsInfo.nlreg(nlregObj0, m1=m1.00, m2=m2.00, 
                                v1=v1.00, v2=v2.00)
  	    Sh <- Shat.nlreg(nlregObj, nlregObj0, m1.1=m1.1, m1.0=m1.0,
                             v1.1=v1.1, v1.0=v1.0)
  	    Sh.t <- Sh[smaller, smaller, drop=FALSE]
  	    qh <- qhat.nlreg(nlregObj, nlregObj0, m1.1=m1.1, v1.1=v1.1)
            s.Sh <- solve(qr(Sh))
            dimnames(s.Sh) <- dimnames(Sh)
  	    Sinvq <- s.Sh[match(offsetName,dimnames(Sh)[[1]]),] %*% qh
  	    lambda.sk.num <-  jh * det(Sh)* det(Sinvq)
  	    lambda.sk.den <-  ih * (xx <- sqrt(abs(det(jt))))
            q.sk[next.idx] <- sign(r[next.idx]) * 
                                abs(lambda.sk.num / lambda.sk.den)
  	    lambda.sk <- abs(q.sk[next.idx] / r[next.idx])
  	    l.sk[next.idx] <- as.numeric(lambda.sk)
  	    np.sk[next.idx] <- log( abs(det(Sh.t) * jh.t / ih.t / xx ) ) / 
                                    r[next.idx]
  	    inf.sk[next.idx] <- (xx <- log(l.sk[next.idx])/r[next.idx]) - 
                                  np.sk[next.idx]
  	    rs.sk[next.idx] <- as.numeric(r[next.idx] + xx)
 	    names(rs.sk[next.idx]) <- offsetName
            rt <- resid(nlregObj0)    
            wt <- sqrt(nlregObj0$weights)      
            Dvt <- matrix(0, nrow=npts, ncol=dp) 
            if(nlregObj$ws$xVar)
              Dvt <- v1.0
            else Dvt[,(db+1):dp] <- v1.0
            Dmt <- m1.0
            if( !nlregObj$ws$missingData )
            {
              Dvt <- Dvt[nlregObj$data$dupl,]
              Dmt <- Dmt[nlregObj$data$dupl,]
            }	
            Dvt <- Dvt/2/wt   
            lvt <- - matrix( rowSums( matvec(Dmvh, rt/wt) ), ncol=1)
            dimnames(lvt) <- list(names(nlregObj$ws$allPar), "")
            if( nlregObj$ws$homVar )
              lvt[names(rc),] <- 0
            aa <- cbind( vecmat( 1/wt^2, Dmt ) + vecmat( rt/wt^2, Dvt[,1:db] ),
                         vecmat( 2*rt/wt^2, Dvt[,(db+1):dp] ) )
#            aa <- cbind( vecmat( 1/wt^2, (Dmt + Dvt[,1:db]) ),
#                         vecmat( 2*rt/wt^2, Dvt[,(db+1):dp] ) )
            lthvt <- Dmvh %*% aa
            dimnames(lthvt) <- list(names(nlregObj$ws$allPar), 
                                    names(nlregObj$ws$allPar))
            dimnames(aa) <- list(rep("", nrow(aa)), 
                                 names(nlregObj$ws$allPar))
            lchvt <-  Dmvh[smaller,,drop=FALSE] %*% aa[,smaller,drop=FALSE]
            AA <- cbind(lvh-lvt, lthvt[, smaller, drop=FALSE])
  	    lambda.fr.num <-  jh * det(AA)
  	    lambda.fr.den <-  det(lthvh) * (xx <- sqrt(abs(det(jt))))
            q.fr[next.idx] <- sign(r[next.idx]) * 
                                abs(lambda.fr.num / lambda.fr.den)
  	    lambda.fr <- abs(q.fr[next.idx] / r[next.idx])
  	    l.fr[next.idx] <- as.numeric(lambda.fr)
            np.fr[next.idx] <- 
              log( abs(det(lchvt[smaller, smaller, drop=FALSE]) /
                       det(lchvh[smaller, smaller, drop=FALSE]) *
                       jh.t/xx) ) / r[next.idx]
  	    inf.fr[next.idx] <- (xx <- log(l.fr[next.idx])/r[next.idx]) - 
                                  np.fr[next.idx]
            rs.fr[next.idx] <- as.numeric(r[next.idx] + xx)
 	    names(rs.fr[next.idx]) <- offsetName                             
	    mat[next.idx,-1] <- c(nlregObj0$coef, nlregObj0$varPar, 
                                  rs.sk[next.idx], rs.fr[next.idx], 
                                  r[next.idx], w[next.idx])
  	    go.on <- ( min(abs(c(w[next.idx], r[next.idx], 
                                 rs.sk[next.idx], rs.fr[next.idx]))) < 2.4 )
  	    if( any(!is.finite(c(w[next.idx], r[next.idx], 
                                 rs.sk[next.idx], rs.fr[next.idx]))) )
  	    {
              x[next.idx] <- mat[next.idx,1] <- NA
  	      w[next.idx] <- r[next.idx] <- NA
              rs.sk[next.idx] <- rs.fr[next.idx] <- NA
  	      q.sk[next.idx] <- l.sk[next.idx] <- NA
              np.sk[next.idx] <- inf.sk[next.idx] <- NA
  	      q.fr[next.idx] <- l.fr[next.idx] <- NA
              np.fr[next.idx] <- inf.fr[next.idx] <- NA
  	      go.on <- FALSE
  	    }
  	    else
  	    {
    	      next.idx <- next.idx + ifelse(first.time, -1, 1)
  	      counter <- counter + 1
              new.start.0 <- c(nlregObj0$coef,   
                nlregObj0$varPar[names(nlregObj0$varPar)!="logs"])
  	      if(followUp)
  	      {
  	        if( !any(match("all.profiles.nlreg", 
                               unlist(lapply(sys.calls(), 
                                             function(x) x[[1]])), 
                               nomatch=FALSE)) )
  	          cat(paste(m$offset, "=", format(of), "\n"))
  	        else  cat(paste(format(of), "\n"))
              }  
            }	
   	  }
  	  else 
  	  {
  	    mat[next.idx,-1] <- c(nlregObj0$coef, nlregObj0$varPar, 
  	                          r[next.idx], w[next.idx])	
  	    go.on <- ( min(abs(c(w[next.idx], r[next.idx]))) < 2.4 )
  	    if( any(!is.finite(c(w[next.idx], r[next.idx]))) )
      	    {
  	      x[next.idx] <- mat[next.idx,1] <- NA
  	      w[next.idx] <- r[next.idx] <- NA
  	      go.on <- FALSE
  	     }
  	    else
  	    {
  	      next.idx <- next.idx + ifelse(first.time, -1, 1)
  	      counter <- counter + 1
              new.start.0 <- c(nlregObj0$coef,  
                nlregObj0$varPar[names(nlregObj0$varPar)!="logs"])
  	      if(followUp)
  	      {
  	        if( !any(match("all.profiles.nlreg", 
                               unlist(lapply(sys.calls(), 
                                             function(x) x[[1]])), 
                               nomatch=FALSE)) )
  	          cat(paste(m$offset, "=", format(of), "\n"))
  	        else  cat(paste(format(of), "\n"))
  	      }
  	    }	
  	  }
        }
        else
        {
          go.on <- TRUE
          next.idx <- next.idx + ifelse(first.time, -1, 1)
  	  counter <- counter + 1
          if(followUp)
  	  {
  	    if( !any(match("all.profiles.nlreg", 
                           unlist(lapply(sys.calls(), 
                                         function(x) x[[1]])), 
                           nomatch=FALSE)) )
  	      cat(paste(m$offset, "=", format(of), "\n"))
  	    else  cat(paste(format(of), "\n"))
            cat("\'nlreg\' did not converge: offset value not considered for profiling\n")  
          }
          else
            cat(paste("\n\'nlreg\' did not converge: offset value",
                      m$offset, "=", format(of), 
                      "not considered for profiling\n"))
        }
      }
      if(first.time)
      {
  	if( next.idx == 0 )
  	  cat("\nlower range too short\n")
  	go.on <- TRUE ; next.idx <- signif/2 + 1
  	new.start.0 <- new.start
  	first.time <- FALSE ; limit <- signif
      }
      else
        if( next.idx == signif )
  	  cat("\nupper range too short\n")	
    }

    if(!followUp && nlregObj$ws$hoa)  cat("\n")
    x <- x[is.finite(w)] ; w <- w[is.finite(w)] ; r <- r[is.finite(r)] 
    if(wantHoa)
    {
      rs.sk <- rs.sk[is.finite(rs.sk)]  
      q.sk <- q.sk[is.finite(q.sk)] ; l.sk <- l.sk[is.finite(l.sk)]	
      np.sk <- np.sk[is.finite(np.sk)] ; inf.sk <- inf.sk[is.finite(inf.sk)]
      rs.fr <- rs.fr[is.finite(rs.fr)]  
      q.fr <- q.fr[is.finite(q.fr)] ; l.fr <- l.fr[is.finite(l.fr)]	
      np.fr <- np.fr[is.finite(np.fr)] ; inf.fr <- inf.fr[is.finite(inf.fr)]
    }
    mat <- mat[is.finite(mat[,2]),]
    dimnames(mat) <- list(c(), c(offsetName, smaller, 
                     if(wantHoa) c("rs.sk", "rs.fr") else NULL, "r", "w"))
    cond1 <- ( x < (mle-omit*s.e.) )   
    cond2 <- ( x > (mle+omit*s.e.) )
    profile <- list(w    = spline(x[cond1 | cond2], w[cond1 | cond2], 
                                  n=n), 
  		    r    = spline(x[cond1 | cond2], r[cond1 | cond2], 
                                  n=n))
   if(wantHoa)
      profile <- c(profile, 
  	           list( rs.sk   = spline(x[cond1 | cond2], 
                                          rs.sk[cond1 | cond2], n=n),
                         q.sk    = spline(x[cond1 | cond2], 
                                          q.sk[cond1 | cond2], n=n),
   			 np.sk	 = spline(x[cond1 | cond2], 
                                          np.sk[cond1 | cond2], n=n),
   		         inf.sk	 = spline(x[cond1 | cond2],
                                          inf.sk[cond1 | cond2], n=n),
                         rs.fr   = spline(x[cond1 | cond2], 
                                          rs.fr[cond1 | cond2], n=n), 
                         q.fr    = spline(x[cond1 | cond2], 
                                          q.fr[cond1 | cond2], n=n),
   			 np.fr	 = spline(x[cond1 | cond2], 
                                          np.fr[cond1 | cond2], n=n),
   		         inf.fr	 = spline(x[cond1 | cond2],
                                          inf.fr[cond1 | cond2], n=n)
                                                                      ))
     profile <- c(profile, 
                   list( mle    = c(mle, se=s.e.), 
                         para   = mat, 
  		         obsI   = if(all) jh.all else NULL,
  	                 hoa    = wantHoa, 
  		         points = counter, 
                         n      = n, 
                         call   = m ))
    attr(profile, "class") <- "nlreg.profile"
    profile	
  }
}

all.profiles.nlreg <- function(fitted, hoa = TRUE, precision = 6, 
                               signif = 30, n = 50, omit = 0.5, 
                               trace = FALSE, call, ...)
{
  cat("\nlong calculation --- may take a while\n")
  m <- match.call()
  new.call <- m
  new.call[1] <- call("profile.nlreg")
  nlregObj <- fitted
  if( nlregObj$ws$hoa )  
    md <- nlregObj$ws$md
  else  
  {
    cat("\ndifferentiating mean function -- may take a while")	
    md <- Dmean(nlregObj)
  }
  if( nlregObj$ws$hoa )  
    vd <- nlregObj$ws$vd
  else  
  {
    cat("\ndifferentiating variance function -- may take a while")
    vd <- Dvar(nlregObj)
  }
  new.call$md <- md
  new.call$vd <- vd
  new.call$all <- TRUE
  followUp <- trace
  wantHoa <- hoa
  allOffsets <- names(nlregObj$ws$allPar)
  npar <- length(allOffsets)
  prof <- list()
  if(followUp && !nlregObj$ws$hoa)  cat("\n")
  for(of in allOffsets)
  {
    if(followUp)
      cat(paste("\nOffset:", of))
    if(nlregObj$ws$hoa)  cat("\n")
      new.call$offset <- of
    prof.temp <- eval(new.call)
    prof <- c(prof, 
              list(list(x = prof.temp$w$x, 
                        y = cbind(w    = prof.temp$w$y, 
                                  r    = prof.temp$r$y, 
                                  rs.sk   = if(wantHoa) prof.temp$rs.sk$y 
                                           else NULL,
                                  rs.fr = if(wantHoa) prof.temp$rs.fr$y 
                                           else NULL),
                                  para = prof.temp$para,
                                  mle  = prof.temp$mle )))
  }
  names(prof) <- allOffsets
  prof <- c(prof, list(npar=npar, obsI=prof.temp$obsI, 
                       hoa=wantHoa, call=call))
  attr(prof, "class") <- "all.nlreg.profiles"
  prof
}		

plot.nlreg.profile <- function(x = stop("nothing to plot"), 
                               alpha = 0.05, add.leg = FALSE, 
                               stats = c("sk", "fr"), 
                               cex = 0.7, cex.lab = 1, cex.axis = 1, 
                               cex.main = 1, lwd1 = 1, lwd2 = 2, 
                               lty1 = "solid", lty2 = "solid", 
                               cl1 = "blue", cl2 = "red", 
                               col = "black", ylim = c(-3,3), ...)
{
  m <- match.call()
  stats <- match.arg(stats)
  if( !match(stats, c("sk","fr"), nomatch=FALSE) )
    stop("choice not valid for \"stats\": either \"sk\" or \"fr\"")
  if(x$hoa) 
    cat(paste("Higher order method used:", 
              ifelse(stats=="sk", "Skovgaard's", "Fraser's"), 
              "r*\n"))
  names(x)[if(stats == "sk") 3:6 else 7:10] <- c("rs", "q", "np", "inf")
  confLimit <- qnorm(1-alpha/2)
  plot(x$w, type="n", ylim=ylim, cex=cex, cex.lab=cex.lab,
       cex.axis=cex.axis, cex.main=cex.main, 
       xlab = paste("parameter of interest", names(x$mle[1])), 
       ylab = "profiles", ...)
  lines(x$w, lty="dashed", col=col, ...)
  lines(x$r, lwd=lwd1, lty=lty1, col=cl1) 
  if(x$hoa)  
    lines(x$rs, lwd=lwd2, lty=lty2, col=cl2)
  abline(h=confLimit, lty="dotted", col=col)
  abline(h=-confLimit, lty="dotted", col=col)
  if(add.leg)
  {
    cat("\nChoose legend position\n")
    if(x$hoa)
      legend(locator(1), c("Wald", "r", 
                           ifelse(stats=="sk", "r* (Sk)", "r* (Fr)")), 
             lwd=c(par("lwd"), lwd1, lwd2),
             lty=c("dashed", lty1, lty2), col=c(col, cl1, cl2))
    else 
      legend(locator(1), c("Wald", "r"), lwd=c(par("lwd"), lwd1),
             lty=c("dashed", lty1), col=c(col, cl1))
  }
  invisible(x)
}			

plot.all.nlreg.profiles <- function(x = stop("nothing to plot"), 
                                    nframe, alpha = 0.05,  
                                    stats = c("sk", "fr"), cex = 0.7, 
                                    cex.lab = 1, cex.axis = 1, 
                                    cex.main = 1, lwd1 = 1, lwd2 = 2, 
                                    lty1 = "solid", lty2 = "solid", 
                                    cl1 = "blue", cl2 = "red", 
                                    col = "black", ylim = c(-3,3), 
                                    ...)
{
  stats <- match.arg(stats)
  if( !match(stats, c("sk","fr"), nomatch=FALSE) )
    stop("choice not valid for \"stats\": either \"sk\" or \"fr\"")
  if(x$hoa) 
    cat(paste("Higher order method used:", 
              ifelse(stats=="sk", "Skovgaard's", "Fraser's"), 
              "r*\n"))
  if(missing(nframe))  
  {
    pAsk <- par("ask")
    par(ask = TRUE) 
    on.exit( par(ask = pAsk) )
  }
  else  
  {
    pMfrow <- par("mfrow")
    par(mfrow = nframe)
    on.exit( par(mfrow = pMfrow) )
  }
  npar <- x$npar
  plot.prof <- function(x, alpha, stats, cex, cex.lab, cex.axis, cex.main, 
                        lwd1, lwd2, lty1, lty2, cl1, cl2, col, ylim, 
                        hoa, ...) 
  {
    confLimit <- qnorm(1-alpha/2)
    matplot(x$x, x$y, type="n", ylim=ylim, lty=lty1, col=col, cex=cex,
            cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main, 
            xlab=paste("parameter of interest", names(x$mle[1])), 
            ylab="profiles", ...)
    lines(x$x, x$y[,1], lty="dashed", col=col, ...)
    lines(x$x, x$y[,2], lwd=lwd1, lty=lty1, col=cl1)
    if(hoa)
      lines(x$x, x$y[,ifelse(stats=="sk", 3, 4)], lwd=lwd2, lty=lty2, 
            col=cl2)
    abline(h=confLimit, lty="dotted", col=col)
    abline(h=-confLimit, lty="dotted", col=col)
  }
  sapply(x[1:npar], plot.prof, alpha=alpha, stats=stats, cex=cex, 
         cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main, 
         lwd1=lwd1, lwd2=lwd2, lty1=lty1, lty2=lty2, cl1=cl1, cl2=cl2, 
         col=col, ylim=ylim, hoa=x$hoa, ...)
  invisible(x)	
}		

print.nlreg.profile <- function(x = stop("nothing to plot"), 
                               alpha = 0.05, add.leg = FALSE, 
                               stats = c("sk", "fr"), 
                               cex = 0.7, cex.lab = 1, cex.axis = 1, 
                               cex.main = 1, lwd1 = 1, lwd2 = 2, 
                               lty1 = "solid", lty2 = "solid", 
                               cl1 = "blue", cl2 = "red", 
                               col = "black", ylim = c(-3,3), ...)
plot.nlreg.profile(x = x, alpha = alpha, add.leg = add.leg, 
                   stats = stats, cex = cex, 
                   cex.lab = cex.lab, cex.axis = cex.axis, 
                   cex.main = cex.main, lwd1 = lwd1, 
                   lwd2 = lwd2, lty1 = lty1, lty2 = lty2, 
                   cl1 = cl1, cl2 = cl2, col = col, 
                   ylim = ylim, ...) 

print.all.nlreg.profiles <- function(x = stop("nothing to plot"), 
                                    nframe, alpha = 0.05,  
                                    stats = c("sk", "fr"), cex = 0.7, 
                                    cex.lab = 1, cex.axis = 1, 
                                    cex.main = 1, lwd1 = 1, lwd2 = 2, 
                                    lty1 = "solid", lty2 = "solid", 
                                    cl1 = "blue", cl2 = "red", 
                                    col = "black", ylim = c(-3,3), 
                                    ...)
plot.all.nlreg.profiles(x = x, nframe = nframe, alpha = alpha,  
                        stats = stats, cex = cex, 
                        cex.lab = cex.lab, cex.axis = cex.axis, 
                        cex.main = cex.main, lwd1 = lwd1, 
                        lwd2 = lwd2, lty1 = lty1, lty2 = lty2, 
                        cl1 = cl1, cl2 = cl2, col = col, 
                        ylim = ylim, ...)

contour.all.nlreg.profiles <- 
      function(x, offset1, offset2, alpha = c(0.1, 0.05), 
               stats = c("sk", "fr"), ret = FALSE,
               plotit = TRUE, drawlabels = FALSE, lwd1 = 1, lwd2 = 1, 
               lty1 = "solid", lty2 = "solid", cl1 = "blue", 
               cl2 = "red", col = "black", pch1 = 1, pch2 = 16, 
               cex = 0.5, ...)
{
  m <- match.call()
  plotIt <- plotit
  if( !plotIt && !ret )
    stop("There is nothing to be returned ...")
  stats <- match.arg(stats)
  if( !match(stats, c("sk", "fr"), nomatch=FALSE) )
    stop("choice not valid for \"stats\": either \"sk\" or \"fr\"")
  which <- ifelse( stats=="sk", "rs.sk", "rs.fr" )
  profilesObj <- x 
  hoa <- profilesObj$hoa
  if(hoa) 
    cat(paste("Higher order method used:", 
              ifelse(stats=="sk", "Skovgaard's", "Fraser's"), 
              "r*\n"))
  npar <- profilesObj$npar
  notGiven <- ( missing(offset1) && missing(offset2) )
  if( !notGiven )
  {
    if( (missing(offset1) && !missing(offset2)) || 
        (!missing(offset1) && missing(offset2)) )
      stop("specify either no or both interest parameters")
    para1 <- paste(substitute(offset1))
    para2 <- paste(substitute(offset2))
    if( para1 == para2 )
      stop("interest parameters must be different")
    para <- c(para1, para2)
    .switch <- match(para, names(profilesObj)[1:npar])
    if( .switch[2] < .switch[1] )
      para <- c(para2, para1)
    npar <- 2
  }
  else  para <- names(profilesObj)[1:npar]
  if(plotIt)
  {
    pPty <- par("pty")
    par(pty = "s") 
    on.exit( par(pty = pPty) )
    if( notGiven )
    {
      split.screen(figs = c(npar,npar))
      on.exit( close.screen(all.screens = TRUE) )
    }
    else 
    {
      pMfrow <- par("mfrow")
      par(mfrow = c(1,2))
      on.exit( par(mfrow = pMfrow) )
    }
  }
  if(ret) 
    contourObj <- list()
  mle.cov <- solve(qr(profilesObj$obsI))
  dimnames(mle.cov) <- dimnames(profilesObj$obsI)
  x.w0.r <- seq(-3, 3, length=10)
  x.w.r <- as.matrix(expand.grid(x.w0.r, x.w0.r))
  quantChi <- sqrt(qchisq(1-alpha, 2))
  quantNorm <- qnorm(1-alpha/2)
  pos1 <- pos2 <- idx <- 0
  for( para1 in para)
  {
    pos1 <- pos1 + 1	
    obj1 <- profilesObj[[para1]]
    mle1 <- obj1$mle
    x.w1 <- seq(-3*mle1[2], 3*mle1[2], length=10)
    if( notGiven && plotIt )
    {			
      screen(pos1+npar*pos2)
      plot(obj1$x, obj1$y[,"w"], type="n", xlab="", ylab="", 
           main=para1, ylim=c(-3,3), ...)
      lines(obj1$x, obj1$y[,"w"], lty="dashed", col=col)
      lines(obj1$x, obj1$y[,"r"], lwd=lwd1, lty=lty1, col=cl1)
      if(hoa)
	lines(obj1$x, obj1$y[,which], lwd=lwd2, lty=lty2, col=cl2)
      abline(h=quantNorm, lty="dotted", col=col)  
      abline(h=-quantNorm, lty="dotted", col=col)
    }
    if( ret && notGiven )
      contourObj <- c(contourObj, 
                      list(list(profile=list(x=obj1$x, y=obj1$y))))
    obj1 <- obj1$para
    if( pos1 != npar )
      for( para2 in para[(pos1+1):npar] )
      {
        idx <- idx + 1
        obj2 <- profilesObj[[para2]]
        mle2 <- obj2$mle
        x.w2 <- seq(-3*mle2[2], 3*mle2[2], length=10)
        x.w <- as.matrix(expand.grid(x.w1, x.w2))
        a.cov <- mle.cov[c(para1,para2), c(para1,para2)]
        a.corr <- sqrt(diag(a.cov))
	a.corr <- vecmat( 1/a.corr, matvec(a.cov, 1/a.corr) )
        val.w <- matrix(diag(x.w %*% solve(a.cov, t(x.w))), 
                        ncol=10)	
        val.wr <- matrix(diag(x.w.r %*% solve(a.corr, t(x.w.r))), 
                         ncol=10)	
        obj2 <- obj2$para
        para2.r1 <- smooth.spline(obj1[,"r"], obj1[,para2], 
                                  all.knots=TRUE)
        r2.para2 <- smooth.spline(obj2[,para2], obj2[,"r"], 
                                  all.knots=TRUE)
        para2.rs1 <- smooth.spline(obj1[,which], obj1[,para2], 
                                   all.knots=TRUE)
        rs2.para2 <- smooth.spline(obj2[,para2], obj2[,which], 
                                   all.knots=TRUE) 
        tau.r1 <- tau.r2 <- para.r1 <- para.r2 <- list()
        tau.rs1 <- tau.rs2 <- para.rs1 <- para.rs2 <- list()
        avg.r <- avg.rs <- list()
        for(qq in seq(along=quantChi))
        {
          cut.r1 <- c(-1, 1)
          temp <- predict(para2.r1, c(-1,1)*quantChi[qq])$y
          temp <- predict(r2.para2, temp)$y
          cut.r2 <- c(temp/quantChi[qq], -1, 1)
          para1.r2 <- smooth.spline(obj2[,"r"], obj2[,para1], 
                                   all.knots=TRUE)
          r1.para1 <- smooth.spline(obj1[,para1], obj1[,"r"], 
                                    all.knots=TRUE)
          temp <- predict(para1.r2, c(-1,1)*quantChi[qq])$y
          temp <- predict(r1.para1, temp)$y
          cut.r1 <- c(cut.r1, temp/quantChi[qq])
          para1.r1 <- smooth.spline(obj1[,"r"], obj1[,para1], 
                                    all.knots=TRUE)
          para2.r2 <- smooth.spline(obj2[,"r"], obj2[,para2], 
                                    all.knots=TRUE)	
          avg.r <- c(avg.r, list((acos(cut.r1) + acos(cut.r2))/2))
          diff.r <- acos(cut.r1) - acos(cut.r2)
          avg.r[[qq]][(diff.r<0) & !is.na(diff.r)] <- 
                             -avg.r[[qq]][(diff.r<0) & !is.na(diff.r)]
          diff.r <- abs(diff.r)
          if( !any(is.na(avg.r[[qq]])) )
          {
            inter.r <- smooth.spline(avg.r[[qq]], diff.r, 
                                     all.knots=TRUE)
            new.diff.r <- predict(inter.r, pi)$y
            avg.r[[qq]] <- c(-pi, avg.r[[qq]], pi)
            diff.r <- c(new.diff.r, diff.r, new.diff.r)
            inter.r <- spline(avg.r[[qq]], diff.r, 300)
            tau.r1 <- c(tau.r1, list(quantChi[qq] * 
                                       cos(inter.r$x+inter.r$y/2)))
            tau.r2 <- c(tau.r2, list(quantChi[qq] * 
                                       cos(inter.r$x-inter.r$y/2)))
            para.r1 <- c(para.r1, list(predict(para1.r1, 
                                               tau.r1[[qq]])$y))
            para.r2 <- c(para.r2, list(predict(para2.r2, 
                                               tau.r2[[qq]])$y))
          }
          else 
          {
            cut1.tmp <- cut.r1*quantChi[qq]
            cut2.tmp <- cut.r2*quantChi[qq]
            tau.r1 <- c(tau.r1, list(cut1.tmp))
            tau.r2 <- c(tau.r2, list(cut2.tmp))
            temp.1 <- predict(para1.r1, cut1.tmp)$y
            temp.2 <- predict(para2.r2, cut2.tmp)$y
            para.r1 <- c(para.r1, list(temp.1))
            para.r2 <- c(para.r2, list(temp.2))
          }	
          if(hoa)
          {
            cut.rs1 <- c(-1, 1)
            temp <- predict(para2.rs1, c(-1,1)*quantChi[qq])$y
            temp <- predict(rs2.para2, temp)$y
            cut.rs2 <- c(temp/quantChi[qq], -1, 1)
            para1.rs2 <- smooth.spline(obj2[,which], obj2[,para1], 
                                       all.knots=TRUE)
            rs1.para1 <- smooth.spline(obj1[,para1], obj1[,which], 
                                       all.knots=TRUE)
            temp <- predict(para1.rs2, c(-1,1)*quantChi[qq])$y
            temp <- predict(rs1.para1, temp)$y
            cut.rs1 <- c(cut.rs1, temp/quantChi[qq])
            avg.rs <- c(avg.rs, list((acos(cut.rs1) + 
                                      acos(cut.rs2))/2))
            diff.rs <- acos(cut.rs1) - acos(cut.rs2)
            avg.rs[[qq]][(diff.rs<0) & !is.na(diff.rs)] <- 
                          -avg.rs[[qq]][(diff.rs<0) & !is.na(diff.rs)]
            diff.rs <- abs(diff.rs)
            if( !any(is.na(avg.rs[[qq]])) )
            {
              inter.rs <- smooth.spline(avg.rs[[qq]], diff.rs, 
                                        all.knots=TRUE)
              new.diff.rs <- predict(inter.rs, pi)$y
              avg.rs[[qq]] <- c(-pi, avg.rs[[qq]], pi)
              diff.rs <- c(new.diff.rs, diff.rs, new.diff.rs)
              inter.rs <- spline(avg.rs[[qq]], diff.rs, 300)
              tau.rs1 <- c(tau.rs1, list(quantChi[qq] * 
                                        cos(inter.rs$x+inter.rs$y/2)))
              tau.rs2 <- c(tau.rs2, list(quantChi[qq] * 
                                        cos(inter.rs$x-inter.rs$y/2)))
              r1.rs1 <- smooth.spline(obj1[,which], obj1[,"r"], 
                                      all.knots=TRUE)
              r2.rs2 <- smooth.spline(obj2[,which], obj2[,"r"], 
                                      all.knots=TRUE)
              tau.rs1[[qq]] <- predict(r1.rs1, tau.rs1[[qq]])$y
              tau.rs2[[qq]] <- predict(r2.rs2, tau.rs2[[qq]])$y
              para.rs1 <- c(para.rs1, list(predict(para1.r1, 
                                                   tau.rs1[[qq]])$y))
              para.rs2 <- c(para.rs2, list(predict(para2.r2, 
                                                   tau.rs2[[qq]])$y))
            }
            else
            {
              cut1.tmp <- cut.rs1*quantChi[qq]
              cut2.tmp <- cut.rs2*quantChi[qq]
              r1.rs1 <- smooth.spline(obj1[,which], obj1[,"r"], 
                                      all.knots=TRUE)
              r2.rs2 <- smooth.spline(obj2[,which], obj2[,"r"], 
                                      all.knots=TRUE)
              temp.1 <- predict(r1.rs1, cut1.tmp)$y
              temp.2 <- predict(r2.rs2, cut2.tmp)$y
              tau.rs1 <- c(tau.rs1, list(temp.1))
              tau.rs2 <- c(tau.rs2, list(temp.2))
              temp.1 <- predict(para1.r1, temp.1)$y
              temp.2 <- predict(para2.r2, temp.2)$y
              para.rs1 <- c(para.rs1, list(temp.1))
              para.rs2 <- c(para.rs2, list(temp.2))
            }	
          }
        }
        if( notGiven && plotIt )
          screen(pos1+idx+npar*pos2)
        xlim <- c(min(c(unlist(tau.r1), 
                        if(hoa) unlist(tau.rs1) else NULL)) , 
                  max(c(unlist(tau.r1), 
                        if(hoa) unlist(tau.rs1) else NULL)))
        xlim <- xlim + c(-0.1,0.1)*diff(xlim)
        ylim <- c(min(c(unlist(tau.r2), 
                        if(hoa) unlist(tau.rs2) else NULL)) , 
                  max(c(unlist(tau.r2), 
                        if(hoa) unlist(tau.rs2) else NULL)))
        ylim <- ylim + c(-0.1,0.1)*diff(ylim)
        if( plotIt )
        {
          if( notGiven )
            plot(tau.r1[[1]], tau.r2[[1]], xlim=xlim, ylim=ylim, 
                 xlab = "", ylab = "", type="n", ...)
          else  
            plot(tau.r1[[1]], tau.r2[[1]], xlim=xlim, ylim=ylim, 
                 xlab = para1, ylab = para2, type="n", 
                 main="N(0,1) scale", ...)
          for(qq in seq(along=quantChi))
          {
            if( !any(is.na(avg.r[[qq]])) )
              lines(tau.r1[[qq]], tau.r2[[qq]], lwd=lwd1, lty=lty1, 
                    col=cl1)
            else
            {
              warning(paste("could not calculate profile pairs sketches for parameters",
                            para1, "and", para2, "!"))
              points(tau.r1[[qq]], tau.r2[[qq]], pch=pch1, col=cl1, 
                     cex=cex)
            }	
            if( hoa )
            {
              if( !any(is.na(avg.rs[[qq]])) )
                lines(tau.rs1[[qq]], tau.rs2[[qq]], lwd=lwd2, 
                      lty=lty2, col=cl2)
              else
              {
                warning(paste("could not calculate profile pairs sketches for parameters",
                               para1, "and", para2, "!"))
                points(tau.rs1[[qq]], tau.rs2[[qq]], pch=pch2, 
                       col=cl2, cex=cex)
              }
            }
          } 
          contour(x.w0.r, x.w0.r, val.wr, levels=qchisq(1-alpha, 2), 
                  add=TRUE, lty="dashed", col=col, 
                  labels=paste(1-alpha), drawlabels=drawlabels, ...)
        }
        if(ret)
          w.r <- list(x.w0.r=x.w0.r, x.w0.r=x.w0.r, val.wr=val.wr, 
                      levels=qchisq(1-alpha,2))
        temp1 <- predict(para2.r1, seq(xlim[1], xlim[2], length=50))
        temp2 <- predict(r2.para2, temp1$y)
        trace1.r <- list(x=temp1$x, y=temp2$y)
        temp1 <- predict(para1.r2, seq(ylim[1], ylim[2], length=50))
        temp2 <- predict(r1.para1, temp1$y)
        trace2.r <- list(x=temp2$y, y=temp1$x)
        if( plotIt )
        {
          lines(trace1.r, lty="solid", col=col)
          lines(trace2.r, lty="solid", col=col)
        }
        if( notGiven && plotIt ) 
          screen(pos1+npar*(pos2+idx))
        xlim <- c(min(c(unlist(para.r1), 
                        if(hoa) unlist(para.rs1) else NULL)) , 
                  max(c(unlist(para.r1), 
                        if(hoa) unlist(para.rs1) else NULL)))
        xlim <- xlim + c(-0.1,0.1)*diff(xlim)
        ylim <- c(min(c(unlist(para.r2), 
                        if(hoa) unlist(para.rs2) else NULL)) , 
                  max(c(unlist(para.r2), 
                        if(hoa) unlist(para.rs2) else NULL)))
        ylim <- ylim + c(-0.1,0.1)*diff(ylim)
        if( plotIt )
        {
          if( notGiven )
            plot(para.r1[[1]], para.r2[[1]], xlim=xlim, ylim=ylim, 
                 xlab = "", ylab = "", type="n", ...)
          else  
            plot(para.r1[[1]], para.r2[[1]], xlim=xlim, ylim=ylim, 
                 xlab = para1, ylab = para2, type="n", 
                 main = "original scale", ...)
          for(qq in seq(along=quantChi))
          {
            if( !any(is.na(avg.r[[qq]])) )
              lines(para.r1[[qq]], para.r2[[qq]], lwd=lwd1, lty=lty1,
                    col=cl1)
            else
            {
              warning(paste("could not calculate profile pairs sketches for parameters",
                            para1, "and", para2, "!"))
              points(para.r1[[qq]], para.r2[[qq]], pch=pch1, col=cl1, 
                     cex=cex)
            }
            if(hoa)
            {
              if( !any(is.na(avg.rs[[qq]])) )
                lines(para.rs1[[qq]], para.rs2[[qq]], lwd=lwd2, 
                      lty=lty2, col=cl2)
              else  
              {
                warning(paste("could not calculate profile pairs sketches for parameters",
                              para1, "and", para2, "!"))
                points(para.rs1[[qq]], para.rs2[[qq]], pch=pch2, 
                       col=cl2, cex=cex)
              }
            }
          }
          contour(x.w1+mle1[1], x.w2+mle2[1], val.w, 
                  levels=qchisq(1-alpha, 2), add=TRUE, lty="dashed", 
                  col=col, labels=paste(1-alpha), 
                  drawlabels=drawlabels, ...)
        }
        if(ret)
          w.para <- list(x.w1.mle=x.w1+mle1[1], x.w2.mle=x.w2+mle2[1],
                         val.w=val.w, levels=qchisq(1-alpha, 2))
        trace1 <- predict(smooth.spline(obj1[,para1], obj1[,para2],
                                        all.knots=TRUE), 
                          seq(xlim[1], xlim[2], length=50))
        temp <- predict(smooth.spline(obj2[,para2], obj2[,para1], 
                        all.knots=TRUE), seq(ylim[1], ylim[2], 
                        length=50))
        trace2 <- list(x=temp$y, y=temp$x)
        if( plotIt )
        {
          lines(trace1, lty="solid", col=col)
          lines(trace2, lty="solid", col=col)
        }
        if(ret)
          if( notGiven )
            contourObj[[length(contourObj)]] <- 
              c(contourObj[[length(contourObj)]], 
                list(list( r         = list(x=tau.r1, y=tau.r2),
                           rs        = if(hoa) 
                                         list(x=tau.rs1, y=tau.rs2) 
                                       else NULL,
                           w         = w.r, 
                           trace1.r  = trace1.r, 
                           trace2.r  = trace2.r,
                           para.r    = list(x=para.r1, y=para.r2),
                           para.rs   = if(hoa) 
                                         list(x=para.rs1, y=para.rs2) 
                                       else NULL,
                           para.w    = w.para,
                           trace1.par = trace1, 
                           trace2.par = trace2)))
          else
            contourObj <- c(contourObj, 
               list( r          = list(x=tau.r1, y=tau.r2),
                     rs         = if(hoa) 
                                    list(x=tau.rs1, y=tau.rs2) 
                                  else NULL,
                     w          = w.r,
                     trace1.r   = trace1.r, 
                     trace2.r = trace2.r,
                     para.r     = list(x=para.r1, y=para.r2),
                     para.rs    = if(hoa) 
                                    list(x=para.rs1, y=para.rs2) 
                                  else NULL,
                     para.w     = w.para,
                     trace1.par = trace1, 
                     trace2.par = trace2,
                     para       = c(para1, para2) ))		
    }
    if(ret && notGiven)
      names(contourObj[[length(contourObj)]]) <- 
                c("profile", if(pos1!=npar) para[(pos1+1):npar]
                             else NULL)
    pos2 <- pos2 + 1
    idx <- 0
  }
  if(ret)
  {
    if( notGiven ) 
      names(contourObj) <- para
    contourObj <- c(contourObj, list(npar=npar, alpha=alpha, hoa=hoa,
                                     stats=stats, fixed=!notGiven, call=m))
    attr(contourObj, "class") <- "nlreg.contours"
    contourObj
  }
  else  invisible(profilesObj)
}
	
plot.nlreg.contours <- function(x, alpha = c(0.1, 0.05), 
                                drawlabels = FALSE, lwd1 = 1, 
                                lwd2 = 1, lty1 = "solid", 
                                lty2 = "solid", cl1 = "blue", 
                                cl2 = "red", col = "black", pch1 = 1, 
                                pch2 = 16, cex = 0.5, ...)
{
  hoa <- x$hoa
  npar <- x$npar
  fixed <- x$fixed
  alpha <- x$alpha
  stats <- x$stats
  which <- ifelse( stats=="sk", "rs.sk", "rs.fr" )
  quantChi <- sqrt(qchisq(1-alpha, 2))
  quantNorm <- qnorm(1-alpha/2)
  if(hoa) 
    cat(paste("Higher order method used:", 
              ifelse(stats=="sk", "Skovgaard's", "Fraser's"), 
              "r*\n"))
  pPty <- par("pty")
  par(pty="s") 
  on.exit( par(pty = pPty) )
  if( !fixed )
  {
    split.screen(figs = c(npar,npar))
    on.exit( close.screen(all.screens = TRUE) )
    para <- names(x)[1:npar]
    pos1 <- pos2 <- idx <- 0
    for( para1 in para )
    {
      pos1 <- pos1 + 1
      obj1 <- x[[para1]]
      screen(pos1+npar*pos2)
      plot(obj1$profile$x, obj1$profile$y[,"w"], type="n", xlab = "", 
           ylab = "", main = para1, ylim=c(-3,3), ...)
      lines(obj1$profile$x, obj1$profile$y[,"w"], lty="dashed", 
            col=col)
      lines(obj1$profile$x, obj1$profile$y[,"r"], lwd=lwd1, lty=lty1,
            col=cl1)
      if(hoa)
	lines(obj1$profile$x, obj1$profile$y[,which], lwd=lwd2, 
              lty=lty2, col=cl2)
      abline(h=quantNorm, lty="dotted", col=col)  
      abline(h=-quantNorm, lty="dotted", col=col)
      if( pos1 != npar )
        for( para2 in para[(pos1+1):npar] )
        {
          idx <- idx + 1
          obj2 <- obj1[[para2]]
          screen(pos1+idx+npar*pos2)
          xlim <- c(min(c(unlist(obj2$r$x), unlist(obj2$rs$x))) , 
                    max(c(unlist(obj2$r$x), unlist(obj2$rs$x))))
          xlim <- xlim + c(-0.1,0.1)*diff(xlim)
          ylim <- c(min(c(unlist(obj2$r$y), unlist(obj2$rs$y))) , 
                    max(c(unlist(obj2$r$y), unlist(obj2$rs$y))))
          ylim <- ylim + c(-0.1,0.1)*diff(ylim)
          plot(obj2$trace1.r, xlab = "", ylab = "", type="l", 
               xlim=xlim, ylim=ylim, col=col, ...)
          lines(obj2$trace2.r, lty="solid", col=col)
          for(qq in 1:length(obj2$r$x))
          {
            if( length(obj2$r$x[[qq]]) > 4 )	
              lines(obj2$r$x[[qq]], obj2$r$y[[qq]], lwd=lwd1, 
                    lty=lty1, col=cl1)
            else
            {
              warning(paste("could not calculate profile pairs sketches for parameters", 
                            para1, "and", para2, "!"))
              points(obj2$r$x[[qq]], obj2$r$y[[qq]], pch=pch1,
                     col=cl1, cex=cex)
            }
            if( hoa )
              if( length(obj2$rs$x[[qq]]) > 4 )
                lines(obj2$rs$x[[qq]], obj2$rs$y[[qq]], lwd=lwd2, 
                      lty=lty2, col=cl2)
              else
              {
                warning(paste("could not calculate profile pairs sketches for parameters", 
                              para1, "and", para2, "!"))
                points(obj2$rs$x[[qq]], obj2$rs$y[[qq]], pch=pch2,
                       col=cl2, cex=cex)
              }
          }
          contour(obj2$w$x.w0.r, obj2$w$x.w0.r, obj2$w$val.wr, 
                  levels=obj2$w$levels, lty="dashed", col=col,
                  add=TRUE, labels=paste(1-alpha),  
                  drawlabels=drawlabels, ...)
          screen(pos1+npar*(pos2+idx))
          xlim <- c(min(c(unlist(obj2$para.r$x), 
                          unlist(obj2$para.rs$x))) , 
                    max(c(unlist(obj2$para.r$x), 
                          unlist(obj2$para.rs$x))))
          xlim <- xlim + c(-0.1,0.1)*diff(xlim)
          ylim <- c(min(c(unlist(obj2$para.r$y), 
                          unlist(obj2$para.rs$y))) , 
                    max(c(unlist(obj2$para.r$y), 
                          unlist(obj2$para.rs$y))))
          ylim <- ylim + c(-0.1,0.1)*diff(ylim)
          plot(obj2$trace1.par$x, obj2$trace1.par$y, xlab = "", 
               ylab = "", xlim=xlim, ylim=ylim, type="l", col=col, 
               ...)
          lines(obj2$trace2.par$x, obj2$trace2.par$y, lty="solid",
                col=col)
          for(qq in 1:length(obj2$para.r$x))
          {
            if( length(obj2$para.r$x[[qq]]) > 4 )
              lines(obj2$para.r$x[[qq]], obj2$para.r$y[[qq]], 
                    lwd=lwd1, col=cl1, lty=lty1)
            else
            {
              warning(paste("could not calculate profile pairs sketches for parameters",
                            para1, "and", para2, "!"))
              points(obj2$para.r$x[[qq]], obj2$para.r$y[[qq]], 
                     pch=pch1, col=cl1, cex=cex)
            }		
            if( hoa )
              if( length(obj2$para.rs$x[[qq]]) > 4 )
                lines(obj2$para.rs$x[[qq]], obj2$para.rs$y[[qq]], 
                      lwd=lwd2, lty=lty2, col=cl2)
              else
              {
                warning(paste("could not calculate profile pairs sketches for parameters",
                              para1, "and", para2, "!"))
                points(obj2$para.rs$x[[qq]], obj2$para.rs$y[[qq]], 
                       pch=pch2, col=cl2, cex=cex)
              }
          }
          contour(obj2$para.w$x.w1.mle, obj2$para.w$x.w2.mle,
                  obj2$para.w$val.w, levels=obj2$para.w$levels,
                  lty="dashed", col=col, add=TRUE,
                  labels=paste(1-alpha), drawlabels=drawlabels,
                  ...)
        }
        pos2 <- pos2 + 1
        idx <- 0		
    }
  }
  else
  {
    para <- x$para
    pMfrow <- par("mfrow")
    par(mfrow=c(1,2))
    on.exit( par(mfrow = pMfrow) )
    xlim <- c(min(c(unlist(x$r$x), unlist(x$rs$x))) , 
              max(c(unlist(x$r$x), unlist(x$rs$x))))
    xlim <- xlim + c(-0.1,0.1)*diff(xlim)
    ylim <- c(min(c(unlist(x$r$y), unlist(x$rs$y))) , 
              max(c(unlist(x$r$y), unlist(x$rs$y))))
    ylim <- ylim + c(-0.1,0.1)*diff(ylim)
    plot(x$trace1.r, xlab = x$para[1], ylab = x$para[2], 
         main = "N(0,1) scale", type="l", xlim=xlim, ylim=ylim, 
         col=col, ...)
    lines(x$trace2.r, lty="solid", col=col, ...)
    for(qq in 1:length(x$r$x))
    {
      if( length(x$r$x[[qq]]) > 4 )
        lines(x$r$x[[qq]], x$r$y[[qq]], lwd=lwd1, lty=lty1, col=cl1)
      else
      {
        warning(paste("could not calculate profile pairs sketches for parameters",
                      para[1], "and", para[2], "!"))
        points(x$r$x[[qq]], x$r$y[[qq]], pch=pch1, col=cl1, cex=cex)
      }
      if( hoa )
        if( length(x$rs$x[[qq]]) > 4 )
          lines(x$rs$x[[qq]], x$rs$y[[qq]], lwd=lwd2, lty=lty2,
                col=cl2)
        else
        {
          warning(paste("could not calculate profile pairs sketches for parameters",
                        para[1], "and", para[2], "!"))
          points(x$rs$x[[qq]], x$rs$y[[qq]], pch=pch2, col=cl2, 
                 cex=cex)
        }
    }
    contour(x$w$x.w0.r, x$w$x.w0.r, x$w$val.wr, levels=x$w$levels, 
            lty="dashed", col=col, add=TRUE, 
            labels=paste(1-alpha), drawlabels=drawlabels, ...)
    xlim <- c(min(c(unlist(x$para.r$x), unlist(x$para.rs$x))) , 
              max(c(unlist(x$para.r$x), unlist(x$para.rs$x))))
    xlim <- xlim + c(-0.1,0.1)*diff(xlim)
    ylim <- c(min(c(unlist(x$para.r$y), unlist(x$para.rs$y))) , 
              max(c(unlist(x$para.r$y), unlist(x$para.rs$y))))
    ylim <- ylim + c(-0.1,0.1)*diff(ylim)
    plot(x$trace1.par, xlab = x$para[1], ylab = x$para[2], 
         main = "original scale", type="l", xlim=xlim, ylim=ylim, 
         col=col, ...)
    lines(x$trace2.par, lty="solid", col=col)
    for(qq in 1:length(x$para.r$x))
    {
      if( length(x$para.r$x[[qq]]) > 4 )
        lines(x$para.r$x[[qq]], x$para.r$y[[qq]], lwd=lwd1, lty=lty1,
              col=cl1)
      else
      {
        warning(paste("could not calculate profile pairs sketches for parameters",
                       para[1], "and", para[2], "!"))
        points(x$para.r$x[[qq]], x$para.r$y[[qq]], pch=pch1, col=cl1,
               cex=cex)
      }	
      if( hoa )
        if( length(x$para.rs$x[[qq]]) > 4 )
          lines(x$para.rs$x[[qq]], x$para.rs$y[[qq]], lwd=lwd2, 
                lty=lty2, col=cl2)
        else
        {
          warning(paste("could not calculate profile pairs sketches for parameters",
                        para[1], "and", para[2], "!"))
          points(x$para.rs$x[[qq]], x$para.rs$y[[qq]], pch=pch2, 
                 col=cl2, cex=cex)
        }
    } 
    contour(x$para.w$x.w1.mle, x$para.w$x.w2.mle, x$para.w$val.w, 
            levels=x$para.w$levels, lty="dashed", col=col,
            add=TRUE, labels=paste(1-alpha), drawlabels=drawlabels,
            ...)
  }
  invisible(x)
}

print.nlreg.contours <- function(x, alpha = c(0.1, 0.05), 
                                drawlabels = FALSE, lwd1 = 1, 
                                lwd2 = 1, lty1 = "solid", 
                                lty2 = "solid", cl1 = "blue", 
                                cl2 = "red", col = "black", pch1 = 1, 
                                pch2 = 16, cex = 0.5, ...)
plot.nlreg.contours(x = x, alpha = alpha, 
                    drawlabels = drawlabels, lwd1 = lwd1, 
                    lwd2 = lwd2, lty1 = lty1, lty2 = lty2, 
                    cl1 = cl1, cl2 = cl2, col = col, 
                    pch1 = pch1, pch2 = pch2, cex = cex, ...)

summary.nlreg.profile <- function(object, alpha = 0.05, twoside = TRUE,
                         digits = NULL, ...)
{
  m <- match.call()
  profile.obj <- object
  quant <- if( twoside )  1-alpha/2
           else  1-alpha
  quant <- qnorm(1-quant) ; quant <- c(quant, -quant)
  wald.ci <- matrix(profile.obj$mle[1] + profile.obj$mle[2]*quant, 
                    ncol=2)
  c1 <- ifelse(wald.ci[,1] < wald.ci[,2], wald.ci[,1], wald.ci[,2])
  c2 <- ifelse(wald.ci[,1] < wald.ci[,2], wald.ci[,2], wald.ci[,1])
  wald.ci <- cbind(c1, c2)
  dimnames(wald.ci) <- list(paste("Wald (", 1-alpha, ")", sep=""), 
                            c("lower", "upper"))
  r.ci <- predict(smooth.spline(profile.obj$r$y, profile.obj$r$x), 
                  quant)
  r.ci <- matrix(r.ci$y, ncol=2)
  c1 <- ifelse(r.ci[,1] < r.ci[,2], r.ci[,1], r.ci[,2])
  c2 <- ifelse(r.ci[,1] < r.ci[,2], r.ci[,2], r.ci[,1])
  r.ci <- cbind(c1, c2)
  dimnames(r.ci) <- list(paste("r (", 1-alpha, ")", sep=""), 
                         c("lower", "upper"))
  if(profile.obj$hoa)
  {	
    rs.sk.ci <- predict(smooth.spline(profile.obj$rs.sk$y, 
                                   profile.obj$rs.sk$x), quant)
    rs.sk.ci <- matrix(rs.sk.ci$y, ncol=2)
    c1 <- ifelse(rs.sk.ci[,1] < rs.sk.ci[,2], rs.sk.ci[,1], rs.sk.ci[,2])
    c2 <- ifelse(rs.sk.ci[,1] < rs.sk.ci[,2], rs.sk.ci[,2], rs.sk.ci[,1])
    rs.sk.ci <- cbind(c1, c2)
    dimnames(rs.sk.ci) <- list(paste("r* - Sk (", 1-alpha, ")", sep=""), 
                            c("lower", "upper"))
    rs.fr.ci <- predict(smooth.spline(profile.obj$rs.fr$y, 
                                   profile.obj$rs.fr$x), quant)
    rs.fr.ci <- matrix(rs.fr.ci$y, ncol=2)
    c1 <- ifelse(rs.fr.ci[,1] < rs.fr.ci[,2], rs.fr.ci[,1], rs.fr.ci[,2])
    c2 <- ifelse(rs.fr.ci[,1] < rs.fr.ci[,2], rs.fr.ci[,2], rs.fr.ci[,1])
    rs.fr.ci <- cbind(c1, c2)
    dimnames(rs.fr.ci) <- list(paste("r* - Fr (", 1-alpha, ")", sep=""), 
                            c("lower", "upper"))

    mat <- rbind(rs.fr.ci, rs.sk.ci, r.ci, wald.ci)
  }
  else  mat <- rbind(r.ci, wald.ci)
  summary.obj <- list( CI = mat )
  if(profile.obj$hoa)
    summary.obj <- c( summary.obj, 
                      list( inf.sk = max(abs(profile.obj$inf.sk$y)), 
                            np.sk = max(abs(profile.obj$np.sk$y)),
                            inf.fr = max(abs(profile.obj$inf.fr$y)), 
                            np.fr = max(abs(profile.obj$np.fr$y)) ) )
  mle <- profile.obj$mle
  offset <- names(mle)[1]
  mle <- matrix(mle, nrow=1)
  dimnames(mle) <- list(offset, c("Estimate", "Std. Error"))
  summary.obj <- c(summary.obj, 
                   list(mle = mle, offset = offset, alpha = alpha,
                        twoside = twoside, points = profile.obj$points,
                        n = profile.obj$n, hoa = profile.obj$hoa, 
                        digits = digits, call = m ))
  attr(summary.obj, "class") <- "summary.nlreg.profile"
  summary.obj
}	

summary.all.nlreg.profiles <- function(object, alpha = 0.05, 
                                       twoside = TRUE, digits = NULL,
                                       ...)
{
  m <- match.call()
  profile.obj <- object
  quant <- if( twoside )  1-alpha/2
           else  1-alpha
  quant <- qnorm(1-quant) ; quant <- c(quant, -quant)
  npar <- profile.obj$npar
  getCI <- function(x, quant, alpha, hoa, ...)
  {
    wald.ci <- matrix(x$mle[1] + x$mle[2]*quant, ncol=2)
    c1 <- ifelse(wald.ci[,1] < wald.ci[,2], wald.ci[,1], wald.ci[,2])
    c2 <- ifelse(wald.ci[,1] < wald.ci[,2], wald.ci[,2], wald.ci[,1])
    wald.ci <- cbind(c1, c2)
    dimnames(wald.ci) <- list(paste("Wald (", 1-alpha, ")", sep=""), 
                              c("lower", "upper"))
    r.ci <- predict(smooth.spline(x$y[,2], x$x), quant)
    r.ci <- matrix(r.ci$y, ncol=2)
    c1 <- ifelse(r.ci[,1] < r.ci[,2], r.ci[,1], r.ci[,2])
    c2 <- ifelse(r.ci[,1] < r.ci[,2], r.ci[,2], r.ci[,1])
    r.ci <- cbind(c1, c2)
    dimnames(r.ci) <- list(paste("r (", 1-alpha, ")", sep=""), 
                           c("lower", "upper"))
    if(hoa)
    {	
      rs.sk.ci <- predict(smooth.spline(x$y[,3], x$x), quant)
      rs.sk.ci <- matrix(rs.sk.ci$y, ncol=2)
      c1 <- ifelse(rs.sk.ci[,1] < rs.sk.ci[,2], rs.sk.ci[,1], rs.sk.ci[,2])
      c2 <- ifelse(rs.sk.ci[,1] < rs.sk.ci[,2], rs.sk.ci[,2], rs.sk.ci[,1])
      rs.sk.ci <- cbind(c1, c2)
      dimnames(rs.sk.ci) <- list(paste("r* - Sk (", 1-alpha, ")", sep=""), 
                              c("lower", "upper"))
      rs.fr.ci <- predict(smooth.spline(x$y[,4], x$x), quant)
      rs.fr.ci <- matrix(rs.fr.ci$y, ncol=2)
      c1 <- ifelse(rs.fr.ci[,1] < rs.fr.ci[,2], rs.fr.ci[,1], rs.fr.ci[,2])
      c2 <- ifelse(rs.fr.ci[,1] < rs.fr.ci[,2], rs.fr.ci[,2], rs.fr.ci[,1])
      rs.fr.ci <- cbind(c1, c2)
      dimnames(rs.fr.ci) <- list(paste("r* - Fr (", 1-alpha, ")", sep=""), 
                              c("lower", "upper"))
      mat <- rbind(rs.fr.ci, rs.sk.ci, r.ci, wald.ci)
    }
    else  mat <- rbind(r.ci, wald.ci)
    mat
  }
  summary.obj <- lapply(profile.obj[1:npar], getCI, alpha=alpha, 
                        quant=quant, hoa=profile.obj$hoa)
  mle <- matrix(unlist(lapply(profile.obj[1:npar], "[[", "mle")), 
                nrow=2)
  offset <- names(profile.obj)[1:npar]
  dimnames(mle) <- list(c("Estimate", "Std. Error"), offset)
  summary.obj <- c(summary.obj, 
                   list(mle= mle, offset = offset, alpha = alpha,
                        twoside = twoside, 
                        hoa = profile.obj$hoa, digits = digits, 
                        call = m ))
  class(summary.obj) <- "summary.all.nlreg.profiles"
  summary.obj
}		
	
print.summary.nlreg.profile <- function(x, 
                    digits = if(!is.null(x$digits)) x$digits
                             else max(3, getOption("digits")-3), ...)
{
  cat(paste("\n", if(x$twoside) "Two-sided" else "One-sided",
            "confidence intervals for", x$offset, "\n"))
  nalpha <- length(x$alpha) 
  idx <- rep(nalpha*(0:4), nalpha) + rep(1:nalpha, rep(5,nalpha))
  xCI <- rbind(x$CI, matrix(rep(NA, 2*nalpha), ncol=2))
  print(xCI[idx,], digits=digits, na.print="")
  print(x$mle, digits=digits)
  cat("\n")
  cat(paste(x$points, "points calculated exactly\n"))
  cat(paste(x$n, "points used in spline interpolation\n"))
  if(x$hoa)
  {
    cat(paste("\nINF (Sk):", signif(x$inf.sk, digits=digits)))
    cat(paste("\nINF (Fr):", signif(x$inf.fr, digits=digits)))
    cat(paste("\nNP (Sk):", signif(x$np.sk, digits=digits)))
    cat(paste("\nNP (Fr):", signif(x$np.fr, digits=digits), "\n"))
  }
}			
	
print.summary.all.nlreg.profiles <- 
        function(x, digits = if(!is.null(x$digits)) x$digits
                             else max(3, getOption("digits")-3), ...)
{
  par.names <- x$offset
  nalpha <- length(x$alpha)
  idxa <- rep(nalpha*(0:4), nalpha) + rep(1:nalpha, rep(5,nalpha))
  cat(paste("\n", if(x$twoside) "Two-sided" else "One-sided",
            "confidence intervals for \n\n"))
  for(idx in 1:length(par.names))
  {
    cat(paste(par.names[idx], ":", 
              signif(x$mle[1,idx], digits=digits), "(", 
              signif(x$mle[2,idx], digits=digits), ")", "\n"))
    xCI <- rbind(x[[idx]], matrix(rep(NA, 2*nalpha), ncol=2))
    print(xCI[idxa,], digits=digits, na.print="")
  }
  invisible(x)
}			

mpl <- function(fitted, ...)  
  UseMethod("mpl")

mpl.nlreg <- function(fitted, offset = NULL, stats = c("sk", "fr"), 
                      control = list(x.tol = 1e-6, rel.tol=1e-6, 
                                     step.min=1/2048, maxit = 100),
	              trace = FALSE, ...)
{
  m <- match.call()
  nlregObj <- fitted
  if( missing(nlregObj) )
    stop("fitted model object is missing, with no default")
  if( nlregObj$ws$homVar )
    stop("models with constant variance are not considered")
  if( !is.null(nlregObj$offset) )
    stop("offset parameter cannot be fixed")
  if( !is.null(offset) )
    if( length(offset) != length(nlregObj$varPar) )
      stop("variance parameters must either be fixed or free")
  stats <- match.arg(stats)
  if( !match(stats, c("sk","fr"), nomatch=FALSE) )
    stop("choice not valid for \"stats\": either \"sk\" or \"fr\"")
  if( !is.null(offset) )
    offset <- offset[match(names(offset), names(nlregObj$varPar))]
  if( !missing(control) )
  {
    if( is.null(control$x.tol) )
      control$x.tol <-  1e-6
    if( is.null(control$rel.tol) )
      control$rel.tol <-  1e-6
    if( is.null(control$step.min) )
      control$step.min <-  1/2048
    if( is.null(control$maxit) )
      control$maxit <- 100
  }        
  n <- length(nlregObj$fitted)
  md <- if( nlregObj$ws$hoa )  nlregObj$ws$md
        else  
        {
          cat("\ndifferentiating mean function -- may take a while")
          Dmean(nlregObj)
        }
  vd <- if( nlregObj$ws$hoa )  nlregObj$ws$vd
        else  
        {
          cat("\ndifferentiating variance function -- may take a while")
          Dvar(nlregObj)
        }
#  attach(nlregObj$data, warn.conflicts = FALSE)		## 20.05.13	
  tmp <- as.list(nlregObj$ws$allPar)
  tmp <- c(nlregObj$data, tmp)
  formals(md) <- tmp
  temp <- do.call("md", tmp)
  m1 <- attr(temp, "gradient") ; m1[!is.finite(m1)] <- 0
  formals(vd) <- tmp
  temp <- do.call("vd", tmp)
  v1 <- attr(temp, "gradient") ; v1[!is.finite(v1)] <- 0
#  detach(nlregObj$data)		## 20.05.13
  if( stats == "fr" )
  {
    rh <- resid(nlregObj)  
    wh <- sqrt(nlregObj$weights)      
    db <- ncol(m1)
    dp <- ncol(v1) + ifelse(!nlregObj$ws$xVar, db, 0) 
    npts <- nrow(m1)
    Dvh <- matrix(0, nrow=npts, ncol=dp) 
    if(nlregObj$ws$xVar)
      Dvh <- v1
    else Dvh[,(db+1):dp] <- v1
    Dmh <- m1
    if( !nlregObj$ws$missingData )
    {
      Dvh <- Dvh[nlregObj$data$dupl,]
      Dmh <- Dmh[nlregObj$data$dupl,]
    }	
    Dvh <- Dvh/2/wh  
    Dmvh <- t(Dvh*rh) 
    Dmvh[1:db,] <- Dmvh[1:db,] + t(Dmh)
  }
  minRegCoef <- function(regCoef, varPar, md, vd, m1, v1, nlregObj, 
                         stats, lastIter=FALSE, is.varPar=FALSE)
  {
    nlreg.temp <- nlregObj
    new.data <- eval(nlreg.temp$call$data)
    if( is.null(new.data) )
    {
      newFrame <- as.list(regCoef)
      names(newFrame) <- names(nlreg.temp$coef)
    }
    else
    {
      class(new.data) <- "list"
      tmp <- as.list(regCoef)
      names(tmp) <- names(nlreg.temp$coef)
      newFrame <- c(new.data, tmp)
    }
    tmp <- as.list(varPar)
    names(tmp) <- names(nlreg.temp$varPar)
    newFrame <- c(newFrame, tmp)
    meanFun <- nlreg.temp$call$formula
    varFun <- nlreg.temp$varFun
    num <- eval( call("-", meanFun[[2]], meanFun[[3]]), 
                 envir=newFrame )
    num <- num^2
    den <- eval( varFun[[2]], envir=newFrame )
    logLik <- sum(log(den)) + sum(num/den)
    if( is.varPar )
    {
      rc <- names(nlregObj$coef) ; vp <- names(nlregObj$varPar) 
      if( stats == "fr" )
        dimnames(Dmvh) <- list(c(rc, vp), rep("", ncol(Dmvh)))
      nlreg.temp$coef <- regCoef ; names(nlreg.temp$coef) <- rc
      nlreg.temp$varPar <- varPar ; names(nlreg.temp$varPar) <- vp
      nlreg.temp$fitted <- eval( meanFun[[3]], envir=newFrame )
      nlreg.temp$weights <- den
      nlreg.temp$ws$allPar <- c(regCoef, varPar)  
      names(nlreg.temp$ws$allPar) <- c(rc, vp)
#      attach(nlreg.temp$data, warn.conflicts = FALSE)		## 20.05.13
      tmp <- as.list(nlreg.temp$ws$allPar)
      tmp <- c(nlreg.temp$data, tmp)
      formals(md) <- tmp
      temp <- do.call("md", tmp)
      m1.t <- attr(temp, "gradient") ; m1.t[!is.finite(m1.t)] <- 0
      m2.t <- attr(temp, "hessian") ; m2.t[!is.finite(m2.t)] <- 0
      formals(vd) <- tmp
      temp <- do.call("vd", tmp)
      v1.t <- attr(temp, "gradient") ; v1.t[!is.finite(v1.t)] <- 0
      v2.t <- attr(temp, "hessian") ; v2.t[!is.finite(v2.t)] <- 0
#      detach(nlreg.temp$data)		## 20.05.13
      if( stats == "sk" )
      {
        S.hat <- Shat.nlreg(nlregObj, nlreg.temp, m1.1=m1, m1.0=m1.t, 
                            v1.1=v1, v1.0=v1.t)
        S.hat <- S.hat[rc, rc]
      }
      else
      {
        rt <- sqrt(num/den)
        wt <- sqrt(den)      
        Dvt <- matrix(0, nrow=npts, ncol=dp) 
        if(nlregObj$ws$xVar)
          Dvt <- v1.t
        else Dvt[,(db+1):dp] <- v1.t
        Dmt <- m1.t
        if( !nlregObj$ws$missingData )
        {
          Dvt <- Dvt[nlregObj$data$dupl,]
          Dmt <- Dmt[nlregObj$data$dupl,]
        }	
        Dvt <- Dvt/2/wt 
        lvt <- - matrix( rowSums( matvec(Dmvh, rt/wt) ), ncol=1)
        dimnames(lvt) <- list(names(nlregObj$ws$allPar), "") 
        aa <- cbind( vecmat( 1/wt^2, Dmt ) + vecmat( rt/wt^2, Dvt[,1:db] ),
                      vecmat( 2*rt/wt^2, Dvt[,(db+1):dp] ) )
#        aa <- cbind( vecmat( 1/wt^2, (Dmt + Dvt[,1:db]) ),
#                     vecmat( 2*rt/wt^2, Dvt[,(db+1):dp] ) ) 
        lchvt <-  Dmvh[rc,,drop=FALSE] %*% aa[,rc,drop=FALSE]
      }
      i.obs <- obsInfo.nlreg(nlreg.temp, m1=m1.t, m2=m2.t, 
                             v1=v1.t, v2=v2.t)
      if( lastIter )
      {
        s.i.obs <- solve(qr(i.obs))
        dimnames(s.i.obs) <- dimnames(i.obs)
        cov.var <- s.i.obs[vp, vp]
      }
      i.obs <- i.obs[rc, rc]
      if( lastIter )
      {
        cov.coef <- solve(qr(i.obs))
        dimnames(cov.coef) <- list(rc, rc)
      }
      lmp <- logLik - log(det(i.obs)) + 
               2*log(det( if(stats == "sk") S.hat else lchvt ))
#      lmp <- logLik - determinant(i.obs)$mod + 
#               2*determinant( if(stats == "sk") S.hat else lchvt )$mod
      if( !lastIter )
	lmp
      else
      {
        mod.tmp <- list(cov.var = cov.var, cov.coef = cov.coef, 
                        lmp = -as.vector(lmp)/2 - n/2*log(2*pi), 
                        lp = -logLik/2 - n/2*log(2*pi) )
        mod.tmp
      }
    }
    else  logLik
  }
  minVarPar <- function(varPar, regCoef, md, vd, m1, v1, nlregObj, stats)
  {
    lmp <- minRegCoef(regCoef=regCoef, varPar=varPar, md=md, vd=vd,
                      m1=m1, v1=v1, nlregObj=nlregObj, stats=stats, 
                      is.varPar=TRUE)
    lmp
  }        
  tol <- control$x.tol ; test1 <- 0 ; test2 <- 0 ; test3 <- 0
  reltol <- control$rel.tol
  maxit <- control$maxit  
  go <- TRUE
  coef.old <- coef.new <- nlregObj$coef
  var.old <- var.new <- if( is.null(offset) )  nlregObj$varPar 
                          else offset
  idx <- 0
  if(trace)  
    if(nlregObj$ws$hoa)  cat("\n")
  else  cat("\n\n")
  while(go) 
  {
    if( is.null(offset) )
    {
      var.new <- optim(par=var.old, fn=minVarPar, regCoef=coef.old, 
                       md=md, vd=vd, m1=m1, v1=v1, nlregObj=nlregObj, 
                       stats=stats, 
#                       control=control, method="BFGS")
                       control=list(maxit=control$maxit, reltol=control$rel.tol), 
                       method="BFGS") 
      var.new <- var.new$p
      test1 <- max(abs((var.old-var.new)/var.new))
      var.old <- var.new
    }
    coef.new <- optim(par=coef.old, fn=minRegCoef, varPar=var.new, 
                      md=md, vd=vd, m1=m1, v1=v1, nlregObj=nlregObj, 
                      stats=stats, 
#                     control=control, method="BFGS")
                      control=list(maxit=control$maxit, reltol=control$rel.tol), 
                      method="BFGS") 
    if( !is.null(offset) )
      lmp <- -coef.new$value/2 - n/2*log(2*pi)
    coef.new <- coef.new$p
    test2 <- max(abs((coef.old-coef.new)/coef.new))
    coef.old <- coef.new
    if( is.null(offset) )
    {
      lmp <- minRegCoef(coef.new, var.new, md, vd, m1, v1, nlregObj, 
                        stats=stats, lastIter=TRUE, is.varPar=TRUE)$lmp
      if( idx==0 )
      {
        lmp.new <- lmp
        test3 <- reltol + 1
      }
      else
      {
        lmp.old <- lmp.new
        lmp.new <- lmp
        test3 <- max(abs((lmp.old-lmp.new)/lmp.new))
      }
    }
    if( is.null(offset) ) idx <- idx + 1
    if( trace && is.null(offset) )
       cat(paste("iteration", idx, 
                 ": modified profile log likelihood =", format(lmp), 
                 "\n"))
    if ( ((test1 < tol) && (test2 < tol)) || (test3 < reltol) ||
         (idx == maxit) )
      go <- FALSE
  }
  if( (maxit>1) && (idx==maxit) )
    warning(paste("\nlinear convergence not obtained in", maxit,
                  "iterations"))
  mod.temp <- minRegCoef(coef.new, var.new, md, vd, m1, v1, nlregObj,
                         stats=stats, lastIter=TRUE, is.varPar=TRUE)
  mplObj <- list( varPar       = if( is.null(offset) ) var.new 
                                   else NULL,
                  coefficients = coef.new,
                  offset       = offset,
                  varParMLE    = nlregObj$varPar,
                  coefMLE      = nlregObj$coef,
                  varParCov    = if( is.null(offset) )  
                                   mod.temp$cov.var else NULL,
                  coefCov      = mod.temp$cov.coef,
                  lmp          = mod.temp$lmp,
                  lp           = mod.temp$lp,
                  stats        = stats, 
                  formula      = nlregObj$call$formula,
                  meanFun      = nlregObj$meanFun,
                  varFun       = nlregObj$varFun,
                  data         = nlregObj$data,
                  nobs         = length(nlregObj$residuals),
                  iter         = if(is.null(offset)) idx else NULL,
                  call         = m,
                  ws           = nlregObj$ws )
  mplObj$ws$iter <- idx
  attr(mplObj, "class") <- c("mpl", "nlreg")
  mplObj
} 

print.mpl <- function(x, digits = max(3, getOption("digits")-3), ...)
{
  stats <- x$stats
  if( !is.null(cl <- x$call) )  
  {
    cat("Formula:\n")
    dput(x$formula)
    cat("Variance function:\n")
    dput(x$varFun)		
  }
  cat(paste("\nHigher order method used:", 
            ifelse(stats=="sk", "Skovgaard's", "Fraser's"), 
            "r*\n"))
  varPar <- x$varPar
  vp <- names(varPar)
  if( !is.null(varPar) )
  {
    varPar <- cbind(varPar, x$varParMLE)
    dimnames(varPar) <- list(vp, c("MMPLE ", "MLE "))		
    cat("\nVariance parameters\n")
    print(varPar, digits=digits, ...)
  }
  else
  {
    cat("\nVariance parameters fixed to:\n")
    print(x$offset, digits=digits, ...)
  }
  coef <- x$coefficients
  rc <- names(coef)
  if( !is.null(coef) )
  {
    coef <- cbind(coef, x$coefMLE)
    dimnames(coef) <- list(rc, c("MMPLE ", "MLE "))		
    cat("\nRegression coefficients\n")
    print(coef, digits=digits, ...)
  }
  else cat("\nNo regression coefficient\n")
  cat("\nTotal number of observations:", x$nobs)
  cat("\nTotal number of parameters:", sum(length(rc)+length(vp)))
  cat("\n-2*Log Lmp", format(-2 * x$lmp, digits=digits), "\n")
  if( !is.null(x$iter) )
    cat("\nAlgorithm converged in", x$iter, 
        if(x$iter > 1) "iterations\n" else "iteration\n")
  invisible(x)
}

summary.mpl <- function(object, correlation = FALSE, digits = NULL, 
                        ...)
{
  mplObj <- object
  mplObj$ws$corr <- correlation
  mplObj$digits <- digits
  class(mplObj) <- "summary.mpl"
  mplObj
}

var2cor.mpl <- function(object, ...)
{
  object <- summary(object)
  var2cor(object)
}

var2cor.summary.mpl <- function(object, ...)
{
  if( !is.null(object$varPar) )
  {
    s.e. <- sqrt(diag( cov <- object$varParCov ))
    corrVarPar <- vecmat( 1/s.e., matvec(cov, 1/s.e.) )
    dimnames(corrVarPar) <- dimnames(object$varParCov)
  }
  else corrVarPar <- NULL
  if( !is.null(object$coefficients) ) 
  {
    s.e. <- sqrt(diag( cov <- object$coefCov ))
    corrCoef <- vecmat( 1/s.e., matvec(cov, 1/s.e.) )
    dimnames(corrCoef) <- dimnames(object$coefCov)
  }
  else corrCoef <- NULL
  list( corrVarPar=corrVarPar, corrCoef=corrCoef )
}

print.summary.mpl <- function(x, corr = FALSE, 
                              digits = if(!is.null(x$digits)) x$digits
                                   else max(3, getOption("digits")-3),
                              ...)
{
  stats <- x$stats
  if(!is.null(cl <- x$call)) 
  {
    cat("Formula:\n")
    dput(x$formula)
    cat("Variance function:\n")
    dput(x$varFun)		
  }
  cat(paste("\nHigher order method used:", 
            ifelse(stats=="sk", "Skovgaard's", "Fraser's"), 
            "r*\n"))
  varPar <- x$varPar
  vp <- names(varPar)
  if( !is.null(varPar) )
  {
    var.se <- sqrt(diag(x$varParCov))
    varPar <- cbind(varPar, x$varParMLE, var.se)
    dimnames(varPar) <- list(vp, c("MMPLE", "MLE", "Std. Error"))
    cat("\nVariance parameters\n")
      print(varPar, digits=digits, ...)
  }
  else
  {
    cat("\nVariance parameters fixed to:\n")
    print(x$offset, digits=digits, ...)
  }
  coef <- x$coefficients
  rc <- names(coef)
  if( !is.null(coef) )
  {
    coef.se <- sqrt(diag(x$coefCov))
    coef <- cbind(coef, x$coefMLE, coef.se)
    dimnames(coef) <- list(rc, c("MMPLE", "MLE", "Std. Error"))	
    cat("\nRegression coefficients\n")
    print(coef, digits=digits, ...)
  }
  else cat("\nNo regression coefficient\n")
  cat("\nTotal number of observations:", x$nobs)
  cat("\nTotal number of parameters:", sum(length(rc)+length(vp)))
  cat("\n-2*Log Lmp", format(-2 * x$lmp, digits=digits), "\n")
  if( missing(corr) )
    corr <- x$ws$corr
  if( corr )
  {
    if( !is.null(varPar) )
    {
      cat("\nCorrelation of variance parameters\n")
      p <- length(vp)
      corr <- vecmat( 1/var.se, matvec(x$varParCov, 1/var.se) )
      dimnames(corr) <- dimnames(x$varParCov)
      ll <- lower.tri(corr)
      corr[ll] <- format(round(corr[ll], digits))
      corr[!ll] <- ""
      print(corr[-1, -p, drop=FALSE], quote=FALSE, digits=digits)
    }
    if( !is.null(rc) ) 
    {
      p <- length(rc)
      if(p > 1) 
      {
        cat("\nCorrelation of regression coefficients\n")
        corr <- vecmat( 1/coef.se, matvec(x$coefCov, 1/coef.se) )
        dimnames(corr) <- dimnames(x$coefCov)
        ll <- lower.tri(corr)
        corr[ll] <- format(round(corr[ll], digits))
        corr[!ll] <- ""
        print(corr[-1, -p, drop=FALSE], quote=FALSE, digits=digits)
      }
    }
  }
  invisible(x)
}
	
obsInfo <- function(object, ...)  
  UseMethod("obsInfo")

obsInfo.nlreg <- function(object, par, mu, v, m1 = NULL, m2 = NULL, 
                          v1 = NULL, v2 = NULL, ...)  
{      
  nlregObj <- object
  rc <- nlregObj$coef 
  vp <- nlregObj$varPar
  of <- nlregObj$offset
  t1 <- nlregObj$data$t1
  t2 <- nlregObj$data$t2
  repl <- nlregObj$data$repl
  if( nlregObj$ws$hoa )
  {
    md <- nlregObj$ws$md
    vd <- nlregObj$ws$vd
  }
  if( missing(par) )  par <- nlregObj$ws$allPar
  if( missing(mu) )  mu <- nlregObj$fitted
  if( missing(v) )  v <- nlregObj$weights
  .probl <- ( nlregObj$ws$homVar && !is.null(of) )
  if( .probl )
    .probl <- ( names(of) =="logs" )		
#  attach(nlregObj$data, warn.conflicts = FALSE)		## 20.05.13
#  on.exit( detach(nlregObj$data) )
  if( !nlregObj$ws$missingData )
  {
    mu <- mu[!duplicated(nlregObj$data$dupl)]
    v <- v[!duplicated(nlregObj$data$dupl)]	        
  }
  if( is.null(m1) || is.null(m2) )
  {
    tmp <- as.list(par)
    tmp <- c(nlregObj$data, tmp)
    temp <- if( nlregObj$ws$hoa ) 
            {     	      
       	      formals(md) <- tmp
              do.call("md", tmp)
            }  
            else  
            {
              tmp.fun <- Dmean(nlregObj)
              formals(tmp.fun) <- tmp
              do.call("tmp.fun", tmp)
            }
    m1 <- attr(temp, "gradient") ; m1[!is.finite(m1)] <- 0
    m2 <- attr(temp, "hessian")  ; m2[!is.finite(m2)] <- 0 
  }
  if( !.probl && ( is.null(v1) || is.null(v2) ) )
  {
    tmp <- as.list(par)
    tmp <- c(nlregObj$data, tmp)
    temp <- if( nlregObj$ws$hoa )
            {
       	      formals(vd) <- tmp
              do.call("vd", tmp)
            }  
            else  
            {
              tmp.fun <- Dvar(nlregObj)
              formals(tmp.fun) <- tmp
              do.call("tmp.fun", tmp)
            }
    v1 <- attr(temp, "gradient") ; v1[!is.finite(v1)] <- 0         
    v2 <- attr(temp, "hessian")  ; v2[!is.finite(v2)] <- 0
  }
  if( nlregObj$ws$homVar )
  {
    v <- v[1]
    if( !.probl )
    {
      v1 <- v1[1,1]
      v2 <- v2[1,1,1]
    }
  }	
  si <- t2 + mu^2*repl - 2*t1*mu
  sir <- 2*m1*(mu*repl - t1)
  if( .probl )
  {
    mat <- array(0, dim=rep(length(rc), 2))
    mat <- colSums( ((mu*repl-t1)/v)*m2[,,,drop=FALSE] ) + 
             crossprod( m1, vecmat(repl/v, m1) )
  }		        
  else
  {
    mat <- array(0, dim=rep(length(c(rc,vp)), 2))
    idx1 <- 1:length(nlregObj$coef)
    if( !is.null(vp) ) idx2 <- length(idx1) +  1:length(vp)    
    mat[idx1,idx1] <- colSums( (1/2*repl/v) * 
                               ( if(nlregObj$ws$xVar) 
                                   v2[,idx1,idx1,drop=FALSE] else 0) + 
                                   ( (mu*repl-t1)/v)*m2[,idx1,idx1,
                                                        drop=FALSE] - 
                                     (1/2*si/v^2) * 
                                     ( if(nlregObj$ws$xVar) 
                                         v2[,idx1,idx1,drop=FALSE] 
                                       else 0) ) + 
                      crossprod(repl/v*m1[,idx1,drop=FALSE],
                                m1[,idx1,drop=FALSE]) +
                      if( nlregObj$ws$xVar)
                        ( crossprod(si/v^3*v1[,idx1,drop=FALSE],
                                    v1[,idx1,drop=FALSE]) - 
                          crossprod(1/2*repl/v^2*v1[,idx1,drop=FALSE],
                                    v1[,idx1,drop=FALSE]) - 
                          crossprod(1/2/v^2*v1[,idx1,drop=FALSE],
                                    sir[,idx1,drop=FALSE]) - 
                          crossprod(1/2/v^2*sir[,idx1,drop=FALSE],
                                   v1[,idx1,drop=FALSE])
                        )
  		      else 0
    if( !is.null(vp) )
    {   	
      mat[idx1,idx2] <- if(nlregObj$ws$xVar)
  	                {
  	                  ( colSums( (1/2*repl/v)*v2[,idx1,idx2,
                                                     drop=FALSE] - 
                                     (1/2*si/v^2)*v2[,idx1,idx2,
                                                     drop=FALSE] ) -
                          crossprod(1/2*repl/v^2*v1[,idx1,drop=FALSE],
                                     v1[,idx2,drop=FALSE]) +
                           crossprod(si/v^3*v1[,idx1,drop=FALSE],
                                     v1[,idx2,drop=FALSE]) ) -
                           crossprod(1/2/v^2*sir[,idx1,drop=FALSE],
                                     v1[,idx2,drop=FALSE] ) 
                        }		  
                        else if( nlregObj$ws$hom )
                               - colSums( 1/2/v^2*sir[,idx1,
                                                      drop=FALSE]*v1 )
  			     else 
                               - crossprod(1/2/v^2*sir[,idx1,
                                                       drop=FALSE],
                                           v1[,,drop=FALSE]) ##
      mat[idx2,idx1] <- t(mat[idx1,idx2])
      mat[idx2,idx2] <- if(nlregObj$ws$xVar)
                          ( colSums( (repl/2/v)*v2[,idx2,idx2,
                                                   drop=FALSE] - 
                                     (si/2/v^2)*v2[,idx2,idx2,
                                                   drop=FALSE] ) - 
                            crossprod(repl/2/v^2*v1[,idx2,drop=FALSE],
                                      v1[,idx2,drop=FALSE]) +
                            crossprod(si/v^3*v1[,idx2,drop=FALSE],
                                      v1[,idx2,drop=FALSE])
                          )
                       else if(!nlregObj$ws$homVar)
                              ( colSums( (repl/2/v)*v2[,,,drop=FALSE] - 
                                       (si/2/v^2)*v2[,,,drop=FALSE] ) - 
                                crossprod(repl/2/v^2*v1[,,drop=FALSE],
                                          v1[,,drop=FALSE]) +
                                crossprod(si/v^3*v1[,,drop=FALSE],
                                          v1[,,drop=FALSE])
                              )
  	                    else
                              ( sum((repl/2/v) - (si/2/v^2))*v2 - 
                                sum(repl/2/v^2)*v1^2 +
                                sum(si/v^3)*v1^2
  			      )
    }
  }
  dimnames(mat) <- list(c(names(rc), names(vp)), 
                        c(names(rc), names(vp)))
  mat
}

expInfo <- function(object, ...)  
  UseMethod("expInfo")		

expInfo.nlreg <- function(object, par, mu, v, m1=NULL, v1=NULL, ...)
{                                       
  nlregObj <- object
  rc <- nlregObj$coef 
  vp <- nlregObj$varPar
  of <- nlregObj$offset
  repl <- nlregObj$data$repl
  if( nlregObj$ws$hoa )
  {
    md <- nlregObj$ws$md
    vd <- nlregObj$ws$vd
  }
  if( missing(par) )  par <- nlregObj$ws$allPar
  if( missing(mu) )  mu <- nlregObj$fitted
  if( missing(v) )  v <- nlregObj$weights
  .probl <- ( nlregObj$ws$homVar && !is.null(of) )
  if(.probl)
    .probl <- ( names(of) =="logs" )
#  attach(nlregObj$data, warn.conflicts = FALSE)		## 20.05.13
#  on.exit( detach(nlregObj$data) )
  if( is.null(m1) )
  {
    tmp <- as.list(par)
    tmp <- c(nlregObj$data, tmp)
    temp <- if( nlregObj$ws$hoa )
            { 
        	  formals(md) <- tmp
              do.call("md", tmp)
            }  
            else  
            {
              tmp.fun <- Dmean(nlregObj, hessian=FALSE)
        	  formals(tmp.fun) <- tmp
              do.call("tmp.fun", tmp)
            }
    m1 <- attr(temp, "gradient")
    m1[!is.finite(m1)] <- 0
  }
  if(!.probl && is.null(v1)) 
  {
    tmp <- as.list(par)
    tmp <- c(nlregObj$data, tmp)
    temp <- if( nlregObj$ws$hoa )
            {
        	  formals(vd) <- tmp
              do.call("vd", tmp)
            }  
            else  
            {
              tmp.fun <- Dvar(nlregObj, hessian=FALSE)
              formals(tmp.fun) <- tmp
              do.call("tmp.fun", tmp)
            }
    v1 <- attr(temp, "gradient")
    v1[!is.finite(v1)] <- 0
  }
  der <- theta.deriv(nlregObj, par=par, mu=mu, v=v, m1=m1, v1=v1)
  if( !nlregObj$ws$missingData )
  {
    mu <- mu[!duplicated(nlregObj$data$dupl)]
    v <- v[!duplicated(nlregObj$data$dupl)]	        
  }	
  if( nlregObj$ws$homVar ) v <- v[1]
  if( .probl )
    der2 <- crossprod( der, vecmat(v*repl, der) )
  else
  {
    der2 <- rbind( cbind( diag(v*repl), 
                          if( nlregObj$ws$homVar) 
                            (xx <- as.vector(repl*2*v*mu)) 
                            else diag(repl*2*v*mu) ),
                   cbind( if( nlregObj$ws$homVar) t(xx) 
                            else diag(repl*2*v*mu), 
                          if( nlregObj$ws$homVar ) 
                            sum(repl*(2*v^2+4*mu^2*v)) 
                            else diag(repl*(2*v^2+4*mu^2*v)) ) )
    der2 <- t(der) %*% der2 %*% der
    if( nlregObj$ws$homVar || !nlregObj$ws$xVar )
    {
       der2[1:length(rc),(length(rc)+1):(length(rc)+length(vp))] <- 0
       der2[(length(rc)+1):(length(rc)+length(vp)),1:length(rc)] <- 0
    }
  }
  dimnames(der2) <- list(c(names(rc), names(vp)), 
                         c(names(rc), names(vp)))
  der2
}

theta.deriv <- function(nlregObj, par, mu, v, m1=NULL, v1=NULL)	
{                                       
  rc <- nlregObj$coef 
  vp <- nlregObj$varPar
  of <- nlregObj$offset
  repl <- nlregObj$data$repl
  if( missing(par) )  par <- nlregObj$ws$allPar
  if( missing(mu) )  mu <- nlregObj$fitted
  if( missing(v) )  v <- nlregObj$weights
  .probl <- ( nlregObj$ws$homVar && !is.null(of) )
  if(.probl)
    .probl <- ( names(of) =="logs" )
#  attach(nlregObj$data, warn.conflicts = FALSE)		## 20.05.13
#  on.exit( detach(nlregObj$data) )
  if( !nlregObj$ws$missingData )
  {
    mu <- mu[!duplicated(nlregObj$data$dupl)]
    v <- v[!duplicated(nlregObj$data$dupl)]	        
  }
  if( nlregObj$ws$homVar ) v <- v[1]
  if( is.null(m1) )
  {
    tmp <- as.list(par)
    tmp <- c(nlregObj$data, tmp)
    temp <- Dmean(nlregObj, hessian=FALSE)
    formals(temp) <- tmp
    temp <- do.call("temp", tmp)
    m1 <- attr(temp, "gradient") ; m1[!is.finite(m1)] <- 0
  }
  if( !.probl && is.null(v1) )
  {
    tmp <- as.list(par)
    tmp <- c(nlregObj$data, tmp)
    temp <- Dvar(nlregObj, hessian=FALSE)
    formals(temp) <- tmp
    temp <- do.call(temp, tmp)
    v1 <- attr(temp, "gradient")  ; v1[!is.finite(v1)] <- 0
  }
  if( .probl )
  {
    mat <- array(0, dim=c(length(repl), length(rc)))
    mat <- 1/v*m1[,,drop=FALSE]
  }		        
  else
  {
    mat <- array(0, dim=c(length(repl)+if( nlregObj$ws$homVar ) 
                                         1 else length(repl), 
                          length(c(rc,vp))))
    idx1 <- 1:length(rc)	
    if( !is.null(vp) )
      idx2 <- length(idx1) + 1:length(vp)
    mat[1:length(repl),idx1] <- 1/v*m1[,idx1,drop=FALSE] - 
                                if( nlregObj$ws$xVar ) 
                                  mu/v^2*v1[,idx1,drop=FALSE] else 0
    if( !is.null(vp) )
      mat[,idx2] <- rbind(as.matrix(-mu/v^2 * if( nlregObj$ws$xVar ) 
                                                v1[,idx2,drop=FALSE] 
                                              else v1[,,drop=FALSE]), 
                          1/2/v^2* if( nlregObj$ws$xVar ) 
                          v1[,idx2,drop=FALSE] else v1[,,drop=FALSE])
    if( nlregObj$ws$xVar )
      mat[(length(repl)+1):(2*length(repl)), idx1] <- 
                                         1/2/v^2*v1[,idx1,drop=FALSE]
  } 
  mat
}

Shat.nlreg <- function(nlregObj1, nlregObj0, par.1, par.0, 
                       mu.1, mu.0, v.1, v.0, m1.1=NULL, m1.0=NULL, 
                       v1.1=NULL, v1.0=NULL)			
{
  rc <- nlregObj1$coef 
  vp <- nlregObj1$varPar
  of <- nlregObj0$offset
  repl <- nlregObj1$data$repl
  .probl <- ( nlregObj1$ws$homVar && is.null(nlregObj1$varPar) )
  if( missing(par.1) )  par.1 <- nlregObj1$ws$allPar
  if( missing(par.0) )  par.0 <- nlregObj0$ws$allPar	
  if( missing(mu.1) )  mu.1 <- nlregObj1$fitted
  if( missing(mu.0) )  mu.0 <- nlregObj0$fitted
  if( missing(v.1) )  v.1 <- nlregObj1$weights
  if( missing(v.0) )  v.0 <- nlregObj0$weights
#  attach(nlregObj1$data, warn.conflicts = FALSE)		## 20.05.13
#  on.exit( detach(nlregObj1$data) )
  if( is.null(m1.1) || is.null(m1.0) )
  {
    md <- if( nlregObj1$ws$hoa )  nlregObj1$ws$md
            else Dmean(nlregObj1, hessian=FALSE)
    tmp <- as.list(nlregObj1$ws$allPar)
    tmp <- c(nlregObj1$data, tmp)
    formals(md) <- tmp
    temp <- do.call("md", tmp)
    m1.1 <- attr(temp, "gradient") ; m1.1[!is.finite(m1.1)] <- 0
    tmp <- as.list(nlregObj0$ws$allPar)
    tmp <- c(nlregObj1$data, tmp)
    formals(md) <- tmp
    temp <- do.call("md", tmp)
    m1.0 <- attr(temp, "gradient") ; m1.0[!is.finite(m1.0)] <- 0
  }
  if( is.null(v1.1) || is.null(v1.0) )
  {
    if( !.probl )
    {
      vd <- if( nlregObj1$ws$hoa )  nlregObj1$ws$vd 
              else Dvar(nlregObj1, hessian=FALSE)
      tmp <- as.list(nlregObj1$ws$allPar)
      tmp <- c(nlregObj1$data, tmp)
      formals(vd) <- tmp
      temp <- do.call("vd", tmp)
      v1.1 <- attr(temp, "gradient") ; v1.1[!is.finite(v1.1)] <- 0
      tmp <- as.list(nlregObj0$ws$allPar)
      tmp <- c(nlregObj1$data, tmp)
      formals(temp) <- tmp
      temp <- do.call("vd", tmp)
      v1.0 <- attr(temp, "gradient") ; v1.0[!is.finite(v1.0)] <- 0
    }
    else v1.0 <- v1.1 <- NULL
  }
  der.1 <- theta.deriv(nlregObj1, par=par.1, mu=mu.1, v=v.1, 
                       m1=m1.1, v1=v1.1)
  der.0 <- theta.deriv(nlregObj1, par=par.0, mu=mu.0, v=v.0, 
                       m1=m1.0, v1=v1.0)
  if( !nlregObj1$ws$missingData )
  {
    mu.1 <- mu.1[!duplicated(nlregObj1$data$dupl)]  
    mu.0 <- mu.0[!duplicated(nlregObj1$data$dupl)]
    v.1 <- v.1[!duplicated(nlregObj1$data$dupl)]  
    v.0 <- v.0[!duplicated(nlregObj1$data$dupl)]
  }	
  if( nlregObj1$ws$homVar ) { v.1 <- v.1[1] ; v.0 <- v.0[1] }
  Sigma <- if( .probl )
             diag(v.1*repl)	
           else
             (rbind( cbind( diag(v.1*repl), 
                            if( nlregObj1$ws$homVar) 
                              (xx <- as.vector(repl*2*v.1*mu.1)) 
                              else diag(repl*2*v.1*mu.1) ),
                     cbind( if( nlregObj1$ws$homVar) t(xx) 
                              else diag(repl*2*v.1*mu.1), 
                            if( nlregObj1$ws$homVar ) 
                              sum(repl*(2*v.1^2+4*mu.1^2*v.1)) 
                              else diag(repl*(2*v.1^2+4*mu.1^2*v.1)) )
             ))
  Shat <- t(der.1) %*% Sigma %*% der.0
  if( !.probl )
    if( nlregObj1$ws$homVar || !nlregObj1$ws$xVar)
      Shat[(length(rc)+1):(length(rc)+length(vp)),1:length(rc)] <- 
        ifelse( abs(Shat[(length(rc)+1):(length(rc)+length(vp)),
                         1:length(rc)]) < 1e-10,
                0, Shat[(length(rc)+1):(length(rc)+length(vp)),
                        1:length(rc)] )
  dimnames(Shat) <- list(c(names(rc), names(vp)), 
                         c(names(rc), names(vp)))
  Shat
}

qhat.nlreg <- function(nlregObj1, nlregObj0, par.1, par.0, 
                       mu.1, mu.0, v.1, v.0,  m1.1=NULL, v1.1=NULL)
{                                       
  rc <- nlregObj1$coef 
  vp <- nlregObj1$varPar
  of <- nlregObj0$offset
  repl <- nlregObj1$data$repl
  .probl <- ( nlregObj1$ws$homVar && is.null(nlregObj1$varPar) )
  if( missing(par.1) )  par.1 <- nlregObj1$ws$allPar
  if( missing(par.0) )  par.0 <- nlregObj0$ws$allPar	
  if( missing(mu.1) )  mu.1 <- nlregObj1$fitted
  if( missing(mu.0) )  mu.0 <- nlregObj0$fitted
  if( missing(v.1) )  v.1 <- nlregObj1$weights
  if( missing(v.0) )  v.0 <- nlregObj0$weights
#  attach(nlregObj1$data, warn.conflicts = FALSE)		## 20.05.13
#  on.exit( detach(nlregObj1$data) )
  if( is.null(m1.1) )
  {
    md <- if( nlregObj1$ws$hoa )  nlregObj1$ws$md
            else Dmean(nlregObj1, hessian=FALSE)
    tmp <- as.list(nlregObj1$ws$allPar)
    tmp <- c(nlregObj1$data, tmp)
    formals(md) <- tmp
    temp <- do.call("md", tmp)
    m1.1 <- attr(temp, "gradient") ; m1.1[!is.finite(m1.1)] <- 0
  }
  if( is.null(v1.1) )
  {
    if( !.probl )
    {
      vd <- if( nlregObj1$ws$hoa )  nlregObj1$ws$vd
              else  Dvar(nlregObj1, hessian=FALSE)
      tmp <- as.list(nlregObj1$ws$allPar)
      tmp <- c(nlregObj1$data, tmp)
      formals(vd) <- tmp
      temp <- do.call("vd", tmp)
      v1.1 <- attr(temp, "gradient") ; v1.1[!is.finite(v1.1)] <- 0
    }
    else v1.1 <- NULL
  }
  der.1 <- theta.deriv(nlregObj1, par=par.1, mu=mu.1, v=v.1, 
                       m1=m1.1, v1=v1.1)
  if( !nlregObj1$ws$missingData )
  {
    mu.1 <- mu.1[!duplicated(nlregObj1$data$dupl)] 
    mu.0 <- mu.0[!duplicated(nlregObj1$data$dupl)]
    v.1 <- v.1[!duplicated(nlregObj1$data$dupl)]  
    v.0 <- v.0[!duplicated(nlregObj1$data$dupl)]
  }	
  if( nlregObj1$ws$homVar ) { v.1 <- v.1[1] ; v.0 <- v.0[1] }
  Sigma <- if( .probl )
             diag(v.1*repl) 
           else
             (rbind( cbind( diag(v.1*repl), 
                            if( nlregObj1$ws$homVar) 
                              (xx <- as.vector(repl*2*v.1*mu.1)) 
                              else diag(repl*2*v.1*mu.1)),
                     cbind( if( nlregObj1$ws$homVar) t(xx) 
                              else diag(repl*2*v.1*mu.1), 
                            if( nlregObj1$ws$homVar ) 
                              sum(repl*(2*v.1^2+4*mu.1^2*v.1)) 
                              else diag(repl*(2*v.1^2+4*mu.1^2*v.1)) )
             ))
  thetah <-  mu.1/v.1
  thetat <-  mu.0/v.0
  if( !.probl )
  {
    thetah <- c(thetah, -1/2/v.1)
    thetat <- c(thetat, -1/2/v.0)
  }
  q <- thetah - thetat
  qhat <- t(der.1) %*% Sigma %*% q
  if( !.probl )
    if( !is.null(nlregObj0$offset) )
    {	
      if( nlregObj1$ws$homVar && 
          (names(nlregObj0$offset) == "logs") )
  	qhat[1:length(rc),1] <- 
          ifelse( abs(qhat[1:length(rc),1]) < 1e-10, 
                  0, qhat[1:length(rc),1] )
    }	
  dimnames(qhat) <- list(c(names(rc), names(vp)), NULL)
  qhat
}
