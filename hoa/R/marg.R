## file marg/R/marg.R, v 1.2-2 2014-03-31
##
##  Copyright (C) 2000-2014 Alessandra R. Brazzale 
##
##  This file is part of the "marg" package for R.  This program is 
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
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
##  USA or look up the web page http://www.gnu.org/copyleft/gpl.html.
##
##  Please send any comments, suggestions or errors found to:
##  Alessandra R. Brazzale, Department of Statistics, University of
##  Padova, Via C. Battisti 241/243, 35121 Padova (PD), Italy.
##  Email: alessandra.brazzale@unipd.it  
##  Web: http://www.stat.unipd.it/~brazzale

rsm.distributions <- structure(.Data = 
             list(g0 = function(y,df,...)  (df+1)/2*log(1+y^2/df),
                  g1 = function(y,df,...)  (df+1)*y/(df+y^2),
                  g2 = function(y,df,...)  (df+1)*(df-y^2)/(df+y^2)^2,
                  g0 = function(y,...)     y+2*log(1+exp(-y)),
                  g1 = function(y,...)     (1-exp(-y))/(1+exp(-y)),
                  g2 = function(y,...)     2*exp(-y)/(1+exp(-y))^2,
                  g0 = function(y,...)     exp(y)-y,
                  g1 = function(y,...)     exp(y)-1,
                  g2 = function(y,...)     exp(y),
                  g0 = function(y,...)     exp(-y)+y,
                  g1 = function(y,...)     -exp(-y)+1,
                  g2 = function(y,...)     exp(-y),
                  g0 = function(y,k,...)   ifelse(abs(y) <= k, y^2/2,
                                                  abs(y)*k - k^2/2),
                  g1 = function(y,k,...)   ifelse(abs(y) <= k, y,
                                                  sign(y)*k),
                  g2 = function(y,k,...)   ifelse(abs(y) <= k, 1, 0)),
                               .Dim = c(3,5),
                               .Dimnames = list(c("g0","g1","g2"),
             c("student","logistic","logWeibull","extreme","Huber")))

student <- function(df=stop("Argument \"df\" is missing, with no default"))
{
  df.values <- 1:100
  df.match <- match(df, df.values, nomatch=FALSE)
  if(!df.match)
    stop("Invalid degrees of freedom for t distribution: integer value between 1 and 100")
  make.family.rsm("student", arg=df)
}

logistic <- function()
{
  make.family.rsm("logistic")
}

logWeibull <- function()
{
  make.family.rsm("logWeibull")
}

extreme <- function()
{
  make.family.rsm("extreme")
}

Huber <- function(k = 1.345)
{
  if(k < 0)
    stop("Invalid tuning constant for Huber-type distribution: must be positive")
  if( (k < 1) || (k > 1.5) )
    warning("Tuning constant for Huber-type distribution preferably between 1 and 1.5")
  make.family.rsm("Huber", arg=k)
}

make.family.rsm <- function(name, arg, ...)
{
  if(is.character(name) && 
     charmatch(name, dimnames(rsm.distributions)[[2]], FALSE))
  {
    g0 <- rsm.distributions[["g0",name]]
    g1 <- rsm.distributions[["g1",name]]
    g2 <- rsm.distributions[["g2",name]]
  }
  else
  {
    obj.deriv <- do.call("eval", 
                         list(expr=parse(text=
                              paste(name, ".distributions", sep=""))))
    g0 <- obj.deriv[["g0",name]]
    g1 <- obj.deriv[["g1",name]]
    g2 <- obj.deriv[["g2",name]]
  }
  family <- list(g0=g0, g1=g1, g2=g2, 
                 df=if(charmatch(name, "student", FALSE)) arg,
                 k=if(charmatch(name, "Huber", FALSE)) arg)
  names(family) <- c("g0", "g1", "g2", "df", "k")
  structure(.Data = c(list(family=name), family), 
            class=c("family.rsm","family"))
}

family.rsm <- function(object, ...)
{
  if( length(object$call$family) > 1 )
    eval(object$call$family)
    else
      do.call(deparse(object$call$family, width.cutoff=500), list())
}

print.family.rsm <- function(x, ...)
{
  cat(x$family, "family\n")
  cat("\n g  : ", deparse(x[["g0"]], width.cutoff=500)) 
  cat("\n g' : ", deparse(x[["g1"]], width.cutoff=500))
  cat("\n g'': ", deparse(x[["g2"]], width.cutoff=500), "\n")
  if(charmatch(x$family, "student", FALSE))
    cat("\n df :", x$df, "\n")
  if(charmatch(x$family, "Huber", FALSE))
    cat("\n k  :", x$k, "\n")
}

dHuber <- function(x, k = 1.345)
{
  cnorm <- (2*pnorm(k)-1) + 2*dnorm(k)/k
  ifelse( abs(x) < k, dnorm(x), 
          exp(k^2/2)*dexp(abs(x), k)/(sqrt(2*pi)*k) )/cnorm
}

pHuber <- function(q, k = 1.345)
{
  cnorm <- (2*pnorm(k)-1) + 2*dnorm(k)/k
  x <- abs(q)
  p <- ifelse( x <= k, dnorm(k)/k + pnorm(k) - pnorm(x),
                      exp(k^2/2)*(1-pexp(x,k))/(sqrt(2*pi)*k) )/cnorm 
  ifelse(q<=0, p, 1-p)
}

qHuber <- function(p, k = 1.345)
{
  cnorm <- sqrt(2*pi)*((2*pnorm(k)-1) + 2*dnorm(k)/k)
  x <- pmin(p, 1-p)
  q <- ifelse( x <= sqrt(2*pi)*dnorm(k)/k/cnorm,
               log(k*cnorm*x)/k - k/2,
               qnorm( abs( 1-pnorm(k) + 
                      x*cnorm/sqrt(2*pi) - dnorm(k)/k) ) )
  ifelse(p < 0.5, q, -q)
}

rHuber <- function(n, k = 1.345)
{
  if( is.na(n) )
    return(NA)
  val <- runif(n)
  qHuber(val, k=k)
}

rsm <- function(formula = formula(data), family = gaussian, 
                data = sys.frame(sys.parent()), dispersion = NULL,
                weights = NULL, subset = NULL, na.action = na.fail,
                offset = NULL, method = "rsm.surv", 
                control = glm.control(maxit=100, trace=FALSE),
                model = FALSE, x = FALSE, y = TRUE, contrasts = NULL,
                ...)
{
  call <- match.call()
  org.call <- call
  if(!charmatch(method, c("model.frame","rsm.fit","rsm.surv"), FALSE))
    stop(c("\nUnimplemented method:", method))
  dist <- as.character(call$family)[1]
  user.def <- FALSE
  if(is.na(dist)) 
    dist <- "gaussian"
  else
  {
    if( pmatch(dist, c("gaussian","student","logWeibull","extreme",
                       "logistic","logExponential","logRayleigh",
                       "Huber"), nomatch=FALSE) )
      dist <- match.arg(dist, c("gaussian","student","logWeibull",
                                "extreme","logistic","logExponential",
                                "logRayleigh","Huber"))
    else user.def <- TRUE
  } 
  if(!user.def && (dist == "gaussian"))
  {
    new.call <- call
    new.call[[1]] <- as.name("glm")
    new.call$family <- dist
    if((method=="rsm.fit") || (method=="rsm.surv"))
      new.call$method <- as.character("glm.fit")
    fit <- eval(new.call, sys.parent())
    if(!(method=="model.frame"))
    {
      fit$call <- call
      attr(fit,"class") <- c("lm", "glm", "rsm")
    }
    return(fit)
  }
  if(missing(data))            
    data <- environment(formula)  
  mf <- match.call(expand.dots = FALSE)
  mf$family <- mf$method <- mf$control <- mf$model <- NULL
  mf$dispersion <- mf$x <- mf$y <- mf$contrasts <- NULL
  mf$... <- NULL
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  if(method == "model.frame")
    return(mf)
  if( !any(match("cond.rsm", unlist(lapply(sys.calls(),
                              function(x) x[[1]])), nomatch=FALSE)) )
  {
    if(!missing(family) &&
       !charmatch(dist, c("logWeibull","extreme","logistic","student",
                          "logExponential","logRayleigh","Huber"),
                  FALSE))
    cat(paste("\n Work with user-defined family:", call$family, "\n"))
  }
  if(!missing(dispersion) && is.numeric(dispersion) && 
          !(dispersion > 0))
    stop("\nScale parameter must be positive")
  mt <- attr(mf, "terms")
  na.act <- attr(mf, "na.action") 
  xvars <- as.character(attr(mt, "variables"))[-1]
  if((yvar <- attr(mt, "response")) > 0)
    xvars <- xvars[-yvar]
  xlev <- if(length(xvars) > 0) 
  {
    xlev <- lapply(mf[xvars], levels)
    xlev[!sapply(xlev, is.null)]
  }
  X <- if(!is.empty.model(mt))
         model.matrix(mt, mf, contrasts)
  Y <- model.response(mf, "numeric")
  if(!is.numeric(Y))
     stop("\nResponse must be numeric")
  nobs <- NROW(Y)
  offset <- model.offset(mf) 
  if(!is.null(offset) && length(offset) != nobs)
    stop("Number of offsets is ", length(offset), ", should equal ",
          nobs, " (number of observations)")
  if(is.null(offset))
    offset <- rep(0, nobs)
  w <- model.weights(mf) 
  wzero <- rep(FALSE, nrow(mf))
  if(!length(w))
    w <- rep(1, nrow(mf))
    else if(any(w < 0))
           stop("\nNegative weights not allowed")
           else
           {
             wzero <- (w == 0)
             Y.org <- Y ; X.org <- X ; offset.org <- offset
             Y <- Y*w ; X <- diag(c(w))%*%X ; offset <- w*offset
             if(any(wzero))
             {
               wpos <- !wzero
               fitted <- resid <- q1 <- q2 <- Y.org
               Y <- Y[wpos] ; X <- as.matrix(X[wpos,])  
               offset <- offset[wpos]
               if( !any(match("cond.rsm", unlist(lapply(sys.calls(),
                               function(x) x[[1]])), nomatch=FALSE)) )
               cat(paste("\n", sum(wzero), 
                   if(sum(wzero)==1) "observation" else "observations",
                   "excluded from the model because of zero weights\n"))
             }
           }
  if(is.empty.model(mt)) method <- "rsm.null"  
  else if(!missing(method))
       {
         if(method=="rsm.surv" && dist=="Huber")
           stop("\"rsm.surv\" method not available for Huber-type distribution")
         if( (dist=="student") && (call$family[[2]] < 3) && 
             (method=="rsm.surv") )
           warning("Change method to \"rsm.fit\" for t distribution")
       }
       else
       {
         method <- if(dist=="Huber" || user.def) "rsm.fit"
                     else "rsm.surv"
       }
  if((dist=="student") && (call$family[[2]] < 3))
    method <- "rsm.fit"
  rsm.fitter <- get(method)
  if(is.call(call$family)) 
    call$family[[1]] <- as.name(dist)
  else 
    call$family <- as.name(dist)
  if(charmatch(dist, "logExponential", FALSE))
  {
    if(!is.null(dispersion))
      stop("Invalid value for scale parameter: must be 1 for logExponential error distribution")
    dispersion <- 1
    call$family <- as.name("logWeibull")
  }
  if(charmatch(dist, "logRayleigh", FALSE))
  {
    if(!is.null(dispersion))
      stop("Invalid value for scale parameter: must be 0.5 for logRayleigh error distribution")
    dispersion <- 0.5
    call$family <- as.name("logWeibull")
  }
  .family <- if(is.call(call$family)) eval(call$family)
            else  do.call(deparse(call$family, width.cutoff=500), list())
  offset4fit <- offset
  score.dispersion <- NULL
  fit <- rsm.fitter(X=X, Y=Y, offset=offset4fit, family=.family,
                    dispersion=dispersion, score.dispersion = score.dispersion,
                    maxit=control$maxit, epsilon=control$epsilon, 
                    trace=control$trace, ...)
  if(any(wzero))
  {
    nas <- is.na(fit$coef)
    fitted[wpos] <- fit$fitted.values/w[wpos]
    fitted[wzero] <- X.org[wzero,!nas]%*%as.vector(fit$coef[!nas]) +
                       if(length(offset.org) > 1) 
                         offset.org[wzero] else 0
    fit$fitted.values <- fitted
    resid[wpos] <- fit$resid
    resid[wzero] <- (Y.org[wzero]-fitted[wzero])/fit$dispersion
    fit$residuals <- resid
    q1[wpos] <- fit$q1 ; q2[wpos] <- fit$q2
    q1[wzero] <- .family$g1(resid[wzero], df=.family$df, k=.family$k)
    q2[wzero] <- .family$g2(resid[wzero], df=.family$df, k=.family$k)
    fit$q1 <- q1 ; fit$q2 <- q2
  }
  else
    fit$fitted.values <- fit$fitted.values/w
  fit$weights <- w
  names(fit$fitted.values) <- names(fit$residuals) <- NULL
  names(fit$q1) <- names(fit$q2) <- NULL
  p <- dim(X)[2]
  if(is.null(p)) 
    p <- 0
  rank <- fit$rank
  df.residuals <- length(if(exists("X.org",frame=sys.nframe())) Y.org
                        else Y)-rank-sum(w==0) - ifelse(fit$fixed,0,1)
  asgn <- attr(if(exists("X.org",frame=sys.nframe())) X.org else X,
               "assign")
  if(rank < p)
  {
    nas <- is.na(fit$coef)
    pasgn <- asgn[!nas]
    if(df.residuals > 0)
       fit$assign.residual <- (rank+1):length(Y)
       fit$R.assign <- pasgn
       fit$x.assign <- asgn
  }
  fit <- c(fit, list( assign = asgn,
                      df.residuals = df.residuals,
                      family = .family,
                      user.def = user.def,
                      dist = dist,
                      formula = formula,
                      data = data,
                      terms = mt,
                      contrasts = attr(X, "contrasts"),
                      offset = offset,	## 25.10.13
                      control = control,
                      call = org.call ))
  if(y) fit$y <- if(exists("Y.org",frame=sys.nframe())) Y.org else Y
  names(fit$y) <- NULL
  if(x) fit$X <- if(exists("X.org",frame=sys.nframe())) X.org else X
  if(model) fit$model <- mf
  attr(fit,"class") <- c("rsm","lm") 
  fit
}

rsm.fit <- function(X, Y, offset, family, dispersion, score.dispersion, 
                    maxit, epsilon, trace, ...)
{
  aux.model <- glm.fit(x=X, y=Y, offset=offset, intercept=FALSE)
  attr(aux.model, "class") <- c("glm","lm")
  start <- aux.model$coef
  if(any(nas <- is.na(start)))
  {
    names(nas)<- dimnames(X)[[2]]
    X <- X[,!nas]
    aux.model <- glm.fit(x=X, y=Y, offset=offset, intercept=FALSE)
    attr(aux.model,"class") <- c("glm","lm")
    start <- aux.model$coef
  }
  is.null.disp <- is.null(dispersion)
  huber.disp <- !is.null.disp && !is.numeric(dispersion)
  if(huber.disp && family$family!="Huber")
    stop("Invalid value for scale parameter: must be positive")
  if( is.null.disp )
    dispersion <- sqrt(summary(aux.model)$dispersion)
  if( huber.disp )
    dispersion <- median(abs(resid(aux.model)))/0.6745
  args <- resid(aux.model)/dispersion
  iter <- 1
  error2 <- error3 <- error4 <- 0
  repeat
  {
    if(trace)
       cat("\n iteration", iter, ":")
    if(family$family == "Huber")
    {
      w.h <- family$g1(args, k=family$k)/args
      aux.model <- glm.fit(x=X, y=Y, weights=w.h, offset=offset,
                           family=gaussian(), intercept=FALSE)
      attr(aux.model, "class") <- c("glm","lm")
      new.start <- coef(aux.model)
    }
    else
    {
      w.1 <- family$g1(args, df=family$df)
      w.2 <- family$g2(args, df=family$df)
      if(any( w.2 < 0 ))
      {
        if(trace)
          cat(" (negative iterative weights encountered!)")
        new.start <- start +
          dispersion*solve(qr(t(X)%*%diag(c(w.2))%*%X))%*%t(X)%*%
                                        matrix(c(w.1),ncol=1)
      }
      else
      {
        y.aux <- w.1/w.2
        aux.model <- glm.fit(x=X, y=y.aux, weights=w.2, family=gaussian(),
                             intercept=FALSE)
        attr(aux.model, "class") <- c("glm","lm")
        new.start <- start + dispersion*coef(aux.model)
      }
    }
    error1 <- abs((new.start-start)/start)
    abs.res <- Y-X%*%new.start-if(!is.null(offset)) offset else 0  
    if(is.null.disp)
    {
      aux.dispersion <- optim( par=log(dispersion), fn=rsm.dispersion,
                               abs.res=abs.res, family=family,
                               arg=switch(family$family, 
                                          student=family$df,
                                          Huber=family$k),
                               method="BFGS",
                               control=list(maxit=maxit, reltol=epsilon,
                                            trace=trace) )
      new.dispersion <- exp(aux.dispersion$par)
      score.dispersion <- list(objective=aux.dispersion$objective,
                               grad.norm=aux.dispersion$grad.norm,
                               message=aux.dispersion$message)
      error2 <- abs((new.dispersion-dispersion)/dispersion)
      dispersion <- new.dispersion
    }
    else if(huber.disp)
           dispersion <- median(abs(abs.res))/0.6745
    start <- new.start
    old.args <- args
    args <- abs.res/dispersion
    if(trace)
    {
      loglik <- -length(args)*log(dispersion) -
                   sum(family$g0(args, df=family$df, k=family$k))
      if(family$family == "student")
      { 
        df <- family$df
        nobs <- length(args)
        loglik <- loglik + nobs*log( gamma((df+1)/2)/gamma(1/2)/
                                     gamma(df/2)/sqrt(df) )
      }
      cat(" log likelihood =", signif(loglik,6))
    }
    error3 <- sqrt(sum((args-old.args)^2)/max(1e-20,sum(old.args^2)))
    if((iter == maxit) || (max(error1, error2, error3) < epsilon))
      break
    iter <- iter + 1
  }
  if(trace) cat("\n")
  if( maxit>1 && iter==maxit )
    warning(paste("\nConvergence not obtained in", maxit, 
                  "iterations"))
  coefs <- rep(NA,length(nas))
  coefs[!nas] <- start
  names(coefs) <- names(nas)
  names(dispersion) <- "dispersion"
  fitted <- as.vector(X%*%start+if(!is.null(offset)) offset else 0) 
  residuals <- (Y-fitted)/dispersion
  w.1 <- family$g1(residuals, df=family$df, k=family$k)
  w.2 <- family$g2(residuals, df=family$df, k=family$k)
  rank <- dim(X)[2]
  if( !any(match("cond.rsm", unlist(lapply(sys.calls(), 
                                           function(x) x[[1]])),
                 nomatch=FALSE)) &&
      !any(match("anova.rsm", unlist(lapply(sys.calls(), 
                                           function(x) x[[1]])),
                 nomatch=FALSE)) )
    if( any(w.2<0) )
      warning("negative iterative weights returned")
  Rnames <- dimnames(X)[[2]]
  nn <- is.null(Rnames)
  if(is.null.disp)
  {
    X <- cbind(X,residuals)
    dimnames(X)[[2]] <- c(Rnames,"scale")
  }
  Rnames <- list(dimnames(X)[[2]],dimnames(X)[[2]])
  R <- as.matrix(t(X)%*%diag(c(w.2))%*%X)
  if(is.null.disp)
    R[rank+1,rank+1] <- R[rank+1,rank+1] + length(residuals)
  attributes(R) <- list(dim=dim(R))
  if(!nn) attr(R, "dimnames") <- Rnames  
  if( any(diag(R) < 0) )
    stop("convergengence not obtained!")
  loglik <- -length(residuals)*log(dispersion) -
               sum(family$g0(residuals, df=family$df, k=family$k))  
  if(family$family == "student")
  {
    df <- family$df
    nobs <- length(residuals)
    loglik <- loglik + nobs*log( gamma((df+1)/2)/gamma(1/2)/
                                 gamma(df/2)/sqrt(df) )
  }
  fit <- list(coefficients = coefs,
              dispersion = dispersion,
              fixed = !is.null.disp,
              residuals = residuals,
              fitted.values = fitted,
              loglik = as.numeric(loglik),
              q1 = w.1,
              q2 = w.2,
              rank = rank,
              R = R,
			  score.dispersion = score.dispersion,
              iter = iter )  
  fit
}

rsm.surv <- function(X, Y, offset, family, dispersion, score.dispersion, 
                     maxit, epsilon, trace, ...)
{
  offset4surv <- offset
  aux.model <- glm.fit(x=X, y=Y, offset=offset4surv, intercept=FALSE)
  attr(aux.model, "class") <- c("glm","lm")
  start <- aux.model$coef
  if(any(nas <- is.na(start)))
  {
    names(nas)<- dimnames(X)[[2]]
    X <- X[,!nas]
    aux.model <- glm.fit(x=X, y=Y, offset=offset4surv, intercept=FALSE)
    attr(aux.model,"class") <- c("glm","lm")
    start <- aux.model$coef
  }
  if( is.null.disp <- is.null(dispersion) )
    dispersion <- log(sqrt(summary(aux.model)$dispersion))
  if(!is.null.disp && !is.numeric(dispersion))
    stop("Invalid value for scale parameter: must be positive")
  dist <- survreg.dist <- family$family
  parms <- NULL
  if(dist == "student")
  {
    survreg.dist <- "t"
    parms <- c(df=family$df)
  }
  if(dist == "logWeibull")
    survreg.dist <- "extreme"
  if(dist == "extreme")
  {
    Y <- -Y ; X <- -X
    offset4surv <- - offset4surv
  }
  sd <- survreg.distributions[[survreg.dist]]
  Y <- exp(Y) ; Y <- Surv(Y) ; Y <- cbind( log(Y[,1]), Y[,2] )
  if(is.null(offset))
    offset4surv <- rep(0,dim(X)[1])
  fixed <- list()
  fixed$scale <- if(!is.null.disp) dispersion else 0
  controlvals <- survreg.control(maxit, epsilon)
  init <- start
  if(is.null.disp)
  {
    init <- c(start, dispersion)
    names(init)[length(init)] <- "dispersion"
  }
  fit <- survreg.fit(x=X, y=Y, offset=offset4surv, init=init, 
                     controlvals=controlvals, dist=survreg.dist, 
                     scale=fixed$scale, parms=parms) 
  Y <- Y[,1]
  if(dist == "extreme")
  {
    Y <- -Y ; X <- -X
    offset4surv <- - offset4surv
  }
  coefs <- rep(NA, length(nas))
  start <- if(is.null.disp) fit$coef[-length(fit$coef)] else fit$coef
  coefs[!nas] <- start
  names(coefs) <- names(nas)
  if(is.null.disp)
    dispersion <- exp(fit$coef[length(fit$coef)])
  names(dispersion) <- "dispersion"
  fitted <- as.vector(X%*%start+offset4surv)
  residuals <- (Y-fitted)/dispersion
  w.1 <- family$g1(residuals, df=family$df, k=family$k)
  w.2 <- family$g2(residuals, df=family$df, k=family$k)
  rank <- dim(X)[2]
  Rnames <- dimnames(X)[[2]]
  nn <- is.null(Rnames)
  if( any(diag(fit$var) < 0) )
    stop("convergengence not obtained!")
  if( (vqr <- qr(fit$var))$rank < nrow(fit$var) )
    stop("singular covariance matrix: convergence not obtained!")
  R <- solve(qr(fit$var))
  R <- as.matrix(R)
  if(is.null.disp)
  {
    R[-(rank+1),-(rank+1)] <- R[-(rank+1),-(rank+1)]*dispersion^2
    R[rank+1,-(rank+1)] <- R[-(rank+1),rank+1] <- 
                           R[-(rank+1),rank+1]*dispersion
    Rnames <- c(Rnames,"scale")
  }
  else
    R <- R*dispersion^2
  Rnames <- list(Rnames,Rnames)
  attributes(R) <- list(dim=dim(R))
  if(!nn) attr(R, "dimnames") <- Rnames
  loglik <- as.numeric(fit$loglik[2])
  fit <- list(coefficients = coefs,
              dispersion = dispersion,
              fixed = !is.null.disp,
              residuals = residuals,
              fitted.values = fitted,
              loglik = loglik,
              q1 = w.1,
              q2 = w.2,
              rank = rank,
              R = R,
              score.dispersion = score.dispersion,
              iter = fit$iter )
  fit
}

rsm.null <- function(X=NULL, Y, offset, family, dispersion, score.dispersion,
                     maxit, epsilon, trace, ...)
{
  is.null.disp <- is.null(dispersion)
  huber.disp <- !is.null.disp && !is.numeric(dispersion)
  iter <- NULL
  nobs <- length(Y)
  if(is.null(offset)) 
    offset <- rep(0, length(Y))
  abs.res <- Y-offset
  R <- NULL
  if(is.null.disp)
  { 
    dispersion <- sqrt(sum(abs.res^2)/nobs)
    aux.dispersion <- optim( par=log(dispersion), fn=rsm.dispersion,
                             abs.res=abs.res, family=family,
                             arg=switch(family$family,
                                        student=family$df,
                                        Huber=family$k),
                             method="BFGS",
                             control=list(maxit=maxit, reltol=epsilon,
                                          trace=trace) )
    dispersion <- exp(aux.dispersion$par)
    score.dispersion <- list(objective=aux.dispersion$objective,
                             grad.norm=aux.dispersion$grad.norm,
                             message=aux.dispersion$message)
    iter <- 1
  }
  else
    if(huber.disp)
      dispersion <- median(abs(abs.res))/0.6745
      
  rank <- 0
  coefs <- NULL
  fitted <- offset
  residuals <- (Y-fitted)/dispersion
  w.1 <- family$g1(residuals, df=family$df, k=family$k)
  w.2 <- family$g2(residuals, df=family$df, k=family$k)
  if(is.null.disp)
  {
    R <- as.matrix(sum(residuals^2*w.2) + nobs)
    attr(R, "dimnames") <- list("scale","scale")
  } 
  loglik <- -nobs*log(dispersion) -
               sum(family$g0(residuals, df=family$df, k=family$k))
  if(family$family == "student")
  {
    df <- family$df
    loglik <- loglik + nobs*log( gamma((df+1)/2)/gamma(1/2)/
                                 gamma(df/2)/sqrt(df) )
  }
  fit <- list( coefficients = coefs,
               dispersion = dispersion,
               fixed = !is.null.disp,
               residuals = residuals,
               fitted.values = fitted,
               loglik = loglik,
               q1 = w.1,
               q2 = w.2,
               rank = rank,
               R = R,
               score.dispersion = score.dispersion,
               iter = iter )
  fit
}

rsm.dispersion <- function(log.dispersion, abs.res, family, arg)
{
  n <- length(abs.res)
  q.1 <- family$g1
  ( sum(q.1(abs.res/exp(log.dispersion), df=arg, k=arg)*abs.res)/
        exp(log.dispersion) - n )^2
}

update.rsm <- function(object, formula., ..., evaluate = TRUE)
{
  call <- object$call
  if( is.null(call) )
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if( !missing(formula.) )
    call$formula <- update.formula(formula(object), formula.)
  if (length(extras) > 0) 
  {
    existing <- !is.na(idx.c <- pmatch(names(extras), names(call)))
    idx.e <- seq(along=names(extras))[existing]
    for( a in seq(along=idx.e) ) 
      call[[idx.c[a]]] <- extras[[idx.e[a]]]
    if( any(!existing) ) 
    {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if (evaluate)
      eval(call, parent.frame())
  else call
}

residuals.rsm <- function(object, type = c("deviance", "pearson", 
                            "response", "r.star", "prob", "deletion"),
                          weighting = "observed", ...)
{
  type <- match.arg(type)
  object.diag <- rsm.diag(object, weighting=weighting)
  switch( type, deviance = object.diag$rd,
                pearson  = object.diag$rp,
                response = object.diag$resid,
                deletion = object.diag$rg,
                r.star   = object.diag$rs,
                prob     = object.diag$rcs )
}

logLik.rsm <- function(object, ...)
{
  if( length(list(...)) )
    warning("extra arguments discarded")
  p <- object$rank
  if( !object$fixed ) 
    p <- p + 1
  val <- object$loglik
  r <- object$residuals
  w <- object$weights
  N <- length(r)
  excl <- w == 0
  if( any(excl) )  
  {
    r <- r[!excl]
    N <- length(r)
  }  
  N0 <- N
  attr(val, "nall") <- N0
  attr(val, "nobs") <- N
  attr(val, "df") <- p 
  class(val) <- "logLik"
  val
}

print.rsm <- function(x, digits = max(3, getOption("digits")-3), ...)
{
  if(!is.null(cl <- x$call)) 
  {
    cat("Call:\n")
    dput(cl)
  }
  coef <- x$coef
  if( is.null(coef) || !any(nas <- is.na(coef)) )
    cat("\nCoefficients:\n")
  else  
  {
    if(is.null(names(coef))) 
      names(coef) <- paste("b", 1:length(coef), sep="")
    coef <- coef[!nas]
    cat("\nCoefficients: (", sum(nas),
        " not defined because of singularities)\n", sep="")
  }
  print(coef, digits=digits, ...)
  cat("\nScale parameter: ", format(x$dispersion, digits=digits),
      if(x$fixed) " (fixed)\n" else "\n")
  cat("\nError distribution: ", x$dist, "\n")
  rank <- x$rank
  if(is.null(rank) && !is.null(coef))
    rank <- sum(!nas)
  nobs <- length(x$residuals) - sum(x$weights==0)
  rdf <- x$df.residuals
  if(is.null(rdf))
    rdf <- nobs - rank - ifelse(x$fixed, 0, 1)
  cat("\nDegrees of Freedom:", nobs, "Total;", rdf, "Residual\n")
  cat("-2*Log Likelihood", format(-2*x$loglik, digits=digits), "\n")
  invisible(x)
}

plot.rsm <- function(x, ...)  rsm.diag.plots(rsmfit=x, ...)

summary.rsm <- function(object, correlation=FALSE, digits=NULL, ...)
{
  coef <- object$coef
  disp <- object$dispersion
  fixed <- object$fixed
  resid <- object$residuals
  wt <- object$weights
  nas <- if(!is.null(coef)) is.na(coef)
         else NULL
  n <- length(resid)
  p <- object$rank
  if(is.null(p)) 
    p <- if(!is.null(coef)) sum(!nas)
         else 0
  rdf <- object$df.residuals
  if(is.null(rdf))
    rdf <- n - p - sum(wt==0) - ifelse(fixed, 0, 1)
  R <- object$R
  if(!is.null(R))
  {
    Rnames <- dimnames(R)
    covun <- solve(qr(R))
    dimnames(covun) <- Rnames
    rowlen <- sqrt(diag(covun))
  }
  else covun <- NULL
  if(!is.null(coef))
  {
    cnames <- names(coef[!nas])
    coef <- matrix(rep(coef[!nas], 4), ncol = 4)
    dimnames(coef) <- list(cnames, c("Estimate", "Std. Error", 
                                     "z value", "Pr(>|z|)"))
    coef[, 2] <- rowlen[1:p] %o% disp
    coef[, 3] <- coef[, 1]/coef[, 2]
    coef[, 4] <- 2 * pnorm( - abs(coef[, 3]) )
  }
  if(!fixed)
  {
    disp <- matrix(c(disp,rowlen[p+1]*disp),ncol=2)
    dimnames(disp) <- list("dispersion", c("Estimate", "Std. Error"))
  }
  if(correlation && !is.null(R))
  {
    correl <- covun * outer(1/rowlen, 1/rowlen)
    dimnames(correl) <- Rnames
  }
  else correl <- NULL
  summary <- list( coefficients = coef,
                   dispersion = disp,
                   fixed = fixed,
                   residuals = resid,
                   cov.unscaled = covun,
                   correlation = correl,
                   family = object$family,
                   loglik = object$loglik,
                   terms = object$terms,
                   df = c(p, rdf),
                   iter = object$iter,
                   nas = nas,
                   call = object$call,
                   digits = digits)
  attr(summary,"class") <- c("summary.rsm")
  summary
}

print.summary.rsm <- function(x, 
                              digits = max(3, getOption("digits")-3), 
                              signif.stars = getOption(
                                                "show.signif.stars"),
                              quote = TRUE, ...)
{
  coef <- x$coef
  if(!is.null(coef))
  {
    nas <- x$nas
    p <- sum(!nas)
    correl <- x$correl
    if(any(nas)) 
    {
      nc <- length(nas)
      cnames <- names(nas)
      coef1 <- array(NA, c(nc, 3), list(cnames, dimnames(coef)[[2]]))
      coef1[!nas, ] <- coef
      coef <- coef1
      if(!is.null(correl)) 
      {
        correl1 <- matrix(NA, nc, nc, dimnames=list(cnames, cnames))
        correl1[!nas,!nas] <- correl[1:p,1:p]
        correl <- correl1
      }
    }
  }
  if(missing(digits) && !is.null(x$digits))
     digits <- x$digits
  cat("Call: ")
  dput(x$call)
  if(!is.null(coef))
  {
    if(any(nas))
      cat("\nCoefficients: (", sum(nas),
          " not defined because of singularities)\n", sep = "")
    else cat("\nCoefficients:\n")
    printCoefmat(coef, digits=digits, signif.stars=signif.stars)
  }
  else 
  {
    cat("\nCoefficients:\n")
    print(coef)
  }
  cat(paste("\nScale parameter: "))
  cat(signif(x$dispersion[1], digits=digits), " (",
      if(x$fixed) "fixed" else signif(x$dispersion[2], 
                                      digits=digits), ")\n")
  cat("\nError distribution: ", x$family$family, "\n")
  df <- x$df
  nobs <- df[2] + df[1] + ifelse(x$fixed, 0, 1)
  cat("\nDegrees of Freedom:", nobs, "Total;", x$df[2], "Residual\n")
  cat("-2*Log Likelihood", format(-2*x$loglik), "\n")
  iter <- x$iter
  if(is.null(iter))
    iter <- 0
  cat("\nNumber of Fisher Scoring Iterations:", format(trunc(iter)),
      "\n")
  if(!is.null(coef))
  {
    if(!is.null(correl))
    {
      p <- dim(correl)[2]
      if(p > 1) 
      {
        cat("\nCorrelation of Coefficients:\n")
        ll <- lower.tri(correl)
        correl[ll] <- format(round(correl[ll], digits))
        correl[!ll] <- ""
        print(correl[-1, -p, drop=FALSE], quote=FALSE, digits=digits)
      }
    }
  }
  invisible(x)
}

vcov.rsm <- function(object, correlation=FALSE, ...)
{
  disp <- object$dispersion
  fixed <- object$fixed
  R <- object$R
  if(!is.null(R))
  {
    Rnames <- dimnames(R)
    covun <- solve(qr(R))
    cova <- disp^2*covun
    dimnames(cova) <- Rnames
    rowlen <- sqrt(diag(covun))
  }
  else cova <- NULL
  if(correlation && !is.null(R))
  {
    correl <- covun * outer(1/rowlen, 1/rowlen)
    dimnames(correl) <- Rnames
  }
   if(correlation)  correl  else  cova
}

anova.rsm <- function(object, ... , test = c("Chisq", "none"))
{
  test <- match.arg(test)
  margs <- function(...)  nargs()
  if(margs(...))
     return(anova.rsmlist(list(object, ...), test=test))
  Terms <- object$terms
  term.labels <- attr(Terms, "term.labels")
  object$call$formula <- object$formula <- formula(Terms)
  nt <- length(term.labels)
  m <- model.frame(object)
  family <- family(object)[[1]]
  y <- model.extract(m, "response")
  loglik <- double(nt+1)
  df.res <- loglik
  if(nt)
  {
    loglik[nt+1] <- -2*object$loglik
    df.res[nt+1] <- object$df.residuals
    fit <- object
    for(iterm in seq(from=nt, to=1, by=-1))
    {
      ww <- fit$weights
      argslist <- list(object=fit,
                       formula=eval(parse(text =
                               paste("~ . -", term.labels[iterm]))))
      fit <- do.call("update", argslist)
      loglik[iterm] <- -2*fit$loglik
      df.res[iterm] <- fit$df.residuals
    }
    dev <- c(NA, -diff(loglik))
    df <- c(NA, -diff(df.res))
  }
  else
  {
    loglik[1] <- -2*object$loglik
    df.res[1] <- dim(y)[1] - attr(Terms, "intercept")
    dev <- df <- as.numeric(NA)
  }
  heading <- c("Analysis of Deviance Table\n",
               paste("Error distribution:", family),
               paste("Scale parameter:", ifelse(object$fixed, "fixed",
                     "free (1 df lost for fitting)")),
               paste("Response: ", as.character(formula(object))[2], 
                     "\n", sep = ""), 
                     "Terms added sequentially (first to last)")
  aod <- data.frame(Df=df, Deviance=dev, "Resid. Df"=df.res,
                    "-2*LL"=loglik, row.names=c("NULL", term.labels), 
                    check.names=FALSE)
  if(test != "none")
    aod <- stat.anova(aod, test, scale=1, df.scale=Inf, n=NROW(y))
  structure(aod, heading=heading, class=c("anova", "data.frame"))
}

anova.rsmlist <- function(object, ... , test = c("Chisq", "none"))
{
  test <- match.arg(test)
  rt <- length(object)
  if(rt == 1)
  {
    object <- object[[1]]
    UseMethod("anova")
  }
  forms <- sapply(object, function(x) as.character(formula(x)))
  subs <- as.logical(match(forms[2,], forms[2,1], FALSE))
  if(!all(subs))
    warning("Some fit objects deleted because response differs from the first model")
  if(sum(subs) == 1)
    stop("The first model has a different response from the rest")
  forms <- forms[,subs]
  object <- object[subs]
  dfres <- sapply(object, "[[", "df.residuals")
  m2loglik <- -2*sapply(object, "[[", "loglik")
#  tl <- lapply(object, labels)
  rt <- length(m2loglik)
  effects <- character(rt)
  dm2loglik <- -diff(m2loglik)
  ddf <- -diff(dfres)
  family <- family(object[[1]])[[1]]
  fixed <- object[[1]]$fixed
  heading <- c("Analysis of Deviance Table\n",
               paste("Error distribution:", family),
               paste("Scale parameter:", ifelse(fixed, "fixed",
                     "free (1 df lost for fitting)")),
               paste("Response:", forms[2, 1], "\n"))
  topnote <- paste("Model ", format(1:rt), ": ", forms[3,],
                   sep="", collapse="\n")
  aod <- data.frame("Resid. Df"=dfres, "-2*LL"=m2loglik, 
                    Df=c(NA, abs(ddf)), 
                    Deviance=c(NA, abs(dm2loglik)), check.names=FALSE)
  if(test != "none")
  {
    nobs <- length(object[[1]]$residuals)
    aod <- stat.anova(aod, test, 1, Inf, nobs)
  }
  structure(aod, heading=c(heading, topnote), 
                 class=c("anova", "data.frame"))
}

rsm.diag <- function(rsmfit, weighting = "observed")
{
##
##  Calculate diagnostics for objects of class "rsm".  The diagnostics
##  calculated are various types of residuals as well as the Cook 
##  statistics and the leverages.  Starting point is the IRLS algorithm.
##  Based upon A.J. Canty's "glm.diag" routine.
##
  family <- rsmfit$family
  user.def <- rsmfit$user.def
  f.name <- family[[1]]
  dispersion <- rsmfit$dispersion
  w <- if(is.null(rsmfit$weights)) rep(1, length(rsmfit$residuals))
           else rsmfit$weights
  wzero <- (w == 0)
  resid <- rsmfit$residuals[!wzero]              ## response residuals
  X <- diag(c(w[!wzero])) %*% model.matrix(rsmfit)[!wzero,]
  dev <- 2*(family$g0(resid, df=rsmfit$family$df, k=rsmfit$family$k) -
            family$g0(0, df=rsmfit$family$df, k=rsmfit$family$k))
  q1 <- as.vector(rsmfit$q1[!wzero])
  q2 <- as.vector(rsmfit$q2[!wzero])
  if(missing(weighting))
  {
    if( f.name == "Huber" )
      weighting <- "score"
    if( ((f.name == "student") || user.def) && any(rsmfit$q2<0) )
      weighting <- "score"
  }
  Q2 <- switch(weighting, "observed" = q2 ,
                          "score"    = q1/(resid) ,
                          "deviance" = q1^2/dev ,
                          "max"      = pmax(q2,0) )
  h <- diag( diag(sqrt(Q2)) %*% X %*% solve(qr(t(X)%*%diag(Q2)%*%X))
                            %*% t(X) %*% diag(sqrt(Q2)) )
  p <- rsmfit$rank
  rd <- sign(resid) * sqrt(dev) / sqrt(1 - h)    ## deviance resid.
  rp <- ifelse( Q2==0, NA, q1/sqrt(Q2*(1 - h)) ) ## pearson/working resid.
  rg <- sign(resid) * sqrt((1-h)*rd^2 + h*rp^2)  ## deletion resid.
  rs <- rd - log(rd/rp)/rd                       ## r.star resid.
  if((f.name == "logExponential") || (f.name == "logRayleigh"))
    f.name <- "logWeibull"
  rcs <- switch(f.name, "student"  = pt(resid, df=rsmfit$family$df),
                        "extreme"  = pweibull(exp(resid), shape=1, 
                                                          scale=1),
                        "logistic" = plogis(resid),
                        "Huber"    = pHuber(resid, 
                                            k=rsmfit$family$k) )
  rcs <- qnorm(rcs)                              ## prob. transf. resid.
  cook <- (h * rp^2)/((1 - h) * p)
  list(resid=resid, rd=rd, rp=rp, rg=rg, rs=rs, rcs=rcs, cook=cook, 
       h=h, dispersion=dispersion)
}

rsm.diag.plots <- function(rsmfit, rsmdiag = NULL, weighting = NULL, 
                           which = NULL, subset = NULL, iden = FALSE, 
                           labels = NULL, ret = FALSE, ...)
{
##
##  Diagnostic plots for objects of class "rsm"
##  Based upon A.J. Canty's "glm.diag.plots" routine.
##
  if(is.null(rsmdiag))
  {
    if(is.null(weighting))
    {
      family <- rsmfit$family
      user.def <- rsmfit$user.def
      f.name <- family[[1]]
      weighting <- "observed"
      if( f.name == "Huber" )
        weighting <- "score"
      if( ((f.name == "student") || user.def) && any(rsmfit$q2<0) )
        weighting <- "score"
    }
    rsmdiag <- rsm.diag(rsmfit, weighting = weighting)
  }
  if(is.null(subset))
     subset <- c(1:length(rsmdiag$h))
  else if(is.logical(subset))
         subset <- (1:length(subset))[subset]
  else if(is.numeric(subset) && all(subset < 0))
          subset <- (1:(length(subset) + length(rsmdiag$h)))[subset]
  else if(is.character(subset))
  {
    if(is.null(labels))
      labels <- subset
    subset <- seq(along = subset)
  }
  w <- if(is.null(rsmfit$weights)) rep(1, length(rsmfit$residuals))
       else rsmfit$weights
  wzero <- (w == 0)
  if( rsmfit$family[[1]] == "Huber" )
     cat(paste("\nNOTE: Diagnostic plots for Huber's least",
               " favourable distribution\n",
               "are only meaningful if the errors are _truly_",
               " Huber-type distributed.\n", sep=""))
  choices <- c("All",
 	       "Response residuals against fitted values",
               "Deviance residuals against fitted values",
               "QQ-plot of deviance residuals",
               "Normal QQ-plot of r* residuals",
               "Cook statistic against h/(1-h)",
               "Cook statistic against observation number\n")
  tmenu <- paste("", choices)
  if( is.null(which) )
    pick <- menu(tmenu, 
                 title="\n Make a plot selection (or 0 to exit)\n")
  else if( !match(which, 2:7, nomatch=FALSE) )
         stop("Choice not valid")
       else pick <- which
  if(pick == 0)
    stop(" No graph required ! ")
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  repeat
  {
    switch(pick,
           "1" = { par(pty="s", mfrow=c(1,1))
##
                   close.screen(all.screens = TRUE)
                   split.screen(c(2, 2))
##
                   screen(1)       ##
##  Plot the response residuals against the fitted values
                   x1 <- rsmfit$fitted.values[!wzero]
                   plot(x1, rsmdiag$resid, xlab = "linear predictor",
                        ylab = "response residuals", ...)
##
                   screen(2)       ##
##  Plot the deviance residuals against the fitted values
                   x2 <- rsmfit$fitted.values[!wzero]
                   plot(x2, rsmdiag$rd, xlab = "linear predictor",
                        ylab = "deviance residuals", ...)
##
                   screen(3)       ##
##  Plot a normal QQ-plot of the deviance residuals
                   y3 <- rsmdiag$rd
                   x3 <- qnorm(ppoints(length(y3)))
                   x3 <- x3[rank(y3)]
                   .lim <- c(min(x3, y3), max(x3,y3))
                   plot(x3, y3, 
                        xlab = paste("quantiles of standard normal"),
                        ylab = "ordered deviance residuals",
                        xlim = .lim, ylim = .lim, ...)
                   abline(0, 1, lty = 2)
##
                   screen(4)       ##
##  Plot a normal QQ-plot of the r* residuals
                   y4 <- rsmdiag$rs
                   x4 <- qnorm(ppoints(length(y4)))[rank(y4)]
       	           .lim <- c(min(x4, y4), max(x4,y4))
                   plot(x4, y4, 
                        xlab = paste("quantiles of standard normal"),
                        ylab = "ordered r* residuals",
                        xlim = .lim, ylim = .lim, ...)
                   abline(0, 1, lty = 2)
##
                   xx <- list(x1, x2, x3, x4)
                   yy <- list(rsmdiag$resid, rsmdiag$resid, y3, y4)
                   if(is.null(labels))
                     labels <- names(model.extract(
                                     model.frame(rsmfit), "response"))
##
                   yes <- iden
                   while(yes)
                   {
##  If interaction with the plots is required then ask the user which plot
##  they wish to interact with and then run identify() on that plot.
##  When the user terminates identify(), reprompt until no further interaction
##  is required and the user inputs a 0.
                     cat("****************************************************\n")
                     cat("Please Input a screen number (1,2,3 or 4)\n")
                     cat("0 will terminate the function \n")
                     num <- scan(n = 1)
                     if((length(num) > 0) && 
                        ((num == 1) || (num == 2) || (num == 3) || 
                         (num == 4))) 
                     {
                       cat(paste("Interactive Identification for screen", 
                                 num, "\n"))
                       cat("left button = Identify, center button = Exit\n")
                       screen(num, new = FALSE)
                       identify(xx[[num]], yy[[num]], labels, ...)
                     }
                     else yes <- FALSE
                   }
##
                   close.screen(all.screens=TRUE)
                   par(ask=TRUE)
                   split.screen(figs=c(1,2))
##
                   screen(1)       ##
##  Plot the Cook statistics against h/(1-h) and draw line to highlight
##  possible influential and high leverage points.
                   hh <- rsmdiag$h/(1 - rsmdiag$h)
                   plot(hh, rsmdiag$cook, xlab = "h/(1-h)",
                        ylab = "Cook statistic", ...)
                   rx <- range(hh)
                   ry <- range(rsmdiag$cook)
                   rank.fit <- rsmfit$rank
                   nobs <- rank.fit + rsmfit$df.residuals +
                           ifelse(rsmfit$fixed,0,1)
                   cooky <- 8/(nobs - 2 * rank.fit)
                   hy <- (2 * rank.fit)/(nobs - 2 * rank.fit)
                   if((cooky >= ry[1]) && (cooky <= ry[2]))
                     abline(h = cooky, lty = 2)
                   if((hy >= rx[1]) && (hy <= rx[2]))
                     abline(v = hy, lty = 2)
##
                   screen(2)       ##
##  Plot the Cook statistics against the observation number in the original
##  data set.
                   plot(subset, rsmdiag$cook, xlab = "case",
                        ylab = "Cook statistic", ...)
                   if((cooky >= ry[1]) && (cooky <= ry[2]))
                     abline(h = cooky, lty = 2)
##                  xx <- list(hh, subset)
                   yy <- list(rsmdiag$cook, rsmdiag$cook)
                   if(is.null(labels))
                     labels <- names(model.extract(
                                     model.frame(rsmfit), "response"))
##
                   yes <- iden
                   while(yes)
                   {
##  If interaction with the plots is required then ask the user which plot
##  they wish to interact with and then run identify() on that plot.
##  When the user terminates identify(), reprompt until no further interaction
##  is required and the user inputs a 0.
                     cat("****************************************************\n")
                     cat("Please Input a screen number (1 or 2)\n")
                     cat("0 will terminate the function \n")
                     num <- scan(n = 1)
                     if((length(num) > 0) && ((num == 1) || (num == 2) ))
                     {
                       cat(paste("Interactive Identification for screen",
                           num, "\n"))
                       cat("left button = Identify, center button = Exit\n")
                       screen(num, new = FALSE)
                       identify(xx[[num]], yy[[num]], labels, ...)
                     }
                     else yes <- FALSE
                   }
##
                   close.screen(all.screens = TRUE)
                   par(ask=FALSE)
                 },
           "2" = { par(pty="s", mfrow=c(1,1))
## Plot the response residuals against the fitted values
                   x1 <- rsmfit$fitted.values[!wzero]
                   plot(x1, rsmdiag$resid, xlab = "linear predictor",
                        ylab = "response residuals", ...)
                   xx <- list(x1)
                    yy <- list(rsmdiag$resid)
##
                 },
           "3" = { par(pty="s", mfrow=c(1,1))
##  Plot the deviance residuals against the fitted values
                   x2 <- rsmfit$fitted.values[!wzero]
                   plot(x2, rsmdiag$rd, xlab = "linear predictor",
                        ylab = "deviance residuals", ...)
                   xx <- list(x2)
                   yy <- list(rsmdiag$rd)
                 },
           "4" = { par(pty="s", mfrow=c(1,1))
##  Plot a normal QQ-plot of the deviance residuals
                   y3 <- rsmdiag$rd
                   x3 <- qnorm(ppoints(length(y3)))
                   x3 <- x3[rank(y3)]
	           .lim <- c(min(x3, y3), max(x3,y3))
                   plot(x3, y3, xlab = "quantiles of standard normal",
                        ylab = "ordered deviance residuals",
                        xlim = .lim, ylim = .lim, ...)
                   abline(0, 1, lty = 2)
                   xx <- list(x3)
                   yy <- list(y3)
                 },
           "5" = {
                   par(pty="s", mfrow=c(1,1))
##  Plot a normal QQ-plot of the r* residuals
                   y4 <- rsmdiag$rs
                   x4 <- qnorm(ppoints(length(y4)))[rank(y4)]
		   .lim <- c(min(x4, y4), max(x4,y4))
                   plot(x4, y4, 
                        xlab = paste("quantiles of standard normal"),
                        ylab = "ordered r* residuals",
                        xlim = .lim, ylim = .lim, ...)
                   abline(0, 1, lty = 2)
                   xx <- list(x4)
                   yy <- list(y4)
                 },
           "6" = { par(pty="s", mfrow=c(1,1))
##  Plot the Cook statistics against h/(1-h) and draw line to highlight
##  possible influential and high leverage points.
                   hh <- rsmdiag$h/(1 - rsmdiag$h)
                   plot(hh, rsmdiag$cook, xlab = "h/(1-h)",
                        ylab = "Cook statistic", ...)
                   rx <- range(hh)
                   ry <- range(rsmdiag$cook)
                   rank.fit <- rsmfit$rank
                   nobs <- rank.fit + rsmfit$df.residuals +
                           ifelse(rsmfit$fixed,0,1)
                   cooky <- 8/(nobs - 2 * rank.fit)
                   hy <- (2 * rank.fit)/(nobs - 2 * rank.fit)
                   if((cooky >= ry[1]) && (cooky <= ry[2]))
                           abline(h = cooky, lty = 2)
                   if((hy >= rx[1]) && (hy <= rx[2]))
                           abline(v = hy, lty = 2)
                   xx <- list(hh)
                   yy <- list(rsmdiag$cook)
                 },
           "7" = { par(pty="s", mfrow=c(1,1))
##  Plot the Cook statistics against the observation number in the original
##  data set.
                   plot(subset, rsmdiag$cook, xlab = "case",
                        ylab = "Cook statistic", ...)
                   ry <- range(rsmdiag$cook)
                   rank.fit <- rsmfit$rank
                   nobs <- rank.fit + rsmfit$df.residuals +
                           ifelse(rsmfit$fixed,0,1)
                   cooky <- 8/(nobs - 2 * rank.fit)
                   if((cooky >= ry[1]) && (cooky <= ry[2]))
                            abline(h = cooky, lty = 2)
                   xx <- list(subset)
                   yy <- list(rsmdiag$cook)
                 } )
    if( (!(pick == 1)) )
    {
      if(is.null(labels))
        labels <- names(model.extract(model.frame(rsmfit),"response"))
      yes <- iden
      while(yes)
      {
##  If interaction with the plots is required then ask the user which plot
##  they wish to interact with and then run identify() on that plot.
##  When the user terminates identify(), reprompt until no further interaction
##  is required and the user inputs a 0.
        cat("****************************************************\n")
        cat("Interactive Identification\n")
        cat("left button = Identify, center button = Exit\n")
        identify(xx[[1]], yy[[1]], labels, ...)
        yes <- FALSE
      }
    }
    if( missing(which) )
       pick <- menu(tmenu,
                    title = "\n Make a plot selection (or 0 to exit)\n")
    if( (pick == 0) || !missing(which) )
    {
      invisible(close.screen(all.screens=TRUE))
      break
    }
  }
  if(ret)  rsmdiag  else invisible()
}

cond <- function(object, offset, ...)  UseMethod("cond")

cond.rsm <- function(object, offset, formula = NULL, family = NULL, 
                     dispersion = NULL,data = sys.frame(sys.parent()),
                     pts = 20, n = max(100, 2*pts), tms = 0.6, 
                     from = NULL, to = NULL,  
                     control = glm.control(...), trace = FALSE, ...)
{
  m <- match.call()
  if(missing(offset)) 
    stop("Argument \"offset\" is missing, with no default")
  if(missing(object))
  {
    if(is.null(formula) || is.null(family))
       stop("Model is missing, with no default")
    if(as.character(m$family)[1] == "gaussian")
      stop("\"Gaussian\" error distribution: no need for conditional inference")
    new.call <- m
    new.call$offset <- new.call$pts <- new.call$n <- NULL
    new.call$tms <- new.call$from <- new.call$to <- NULL
    new.call$trace <- NULL
    new.call[[1]] <- as.name("rsm")
    object <- eval(new.call, parent.frame()) 	
  }
  else
  {
    oc <- attr(object, "class")
    if(is.null(oc) || !match("rsm", oc, nomatch=FALSE))
      stop("Invalid argument: not an \"rsm\" object")
    if(object$family[1] == "gaussian")
      stop("\"Gaussian\" error distribution: no need for conditional inference")
  }
  if( object$family[[1]] == "Huber" )
    cat(paste("\nNOTE: Higher order inference for Huber's least",
              " favourable distribution\n",
              "is only meaningful if the errors are _truly_",
              " Huber-type distributed.\n", sep=""))
  mf <- model.frame.default(object, data=object$data)
  Terms <- attr(mf, "terms")
  is.fixed <- object$fixed
  is.empty <- is.empty.model(Terms)
  if(is.empty && is.fixed)
    stop("Invalid model: all parameters are fixed")
  is.scalar <- if(is.empty) !is.fixed 
               else (is.fixed && (dim(model.matrix(object))[2] == 1))
  offsetName <- if(!is.character(m$offset)) deparse(m$offset, 
                                                    width.cutoff=500)
                else m$offset
  .offsetName <- offsetName   
  if((offsetName == "1") && (attr(Terms, "intercept") == 0))
    stop(paste("Invalid argument for \"offset\": intercept not included in original model"))
  if(!match(offsetName, c("scale", "1", if(!is.empty)
                          dimnames(model.matrix(object))[[2]]),  
            nomatch=FALSE))
    stop(paste("Invalid argument for \"offset\": variable not included in original model"))
  if(offsetName == "scale")
  {
    if(is.fixed)
      stop("Invalid argument for \"offset\": scale parameter is already fixed")
  }
  else
  {
    if(is.empty)
      stop("Invalid argument for \"offset\": model is empty")
    switch(offsetName, 
           "1"= { 
                  nobs <- length(model.response(mf))
                  .offset <- rep(1, nobs) 
                }, 
                { 
                  if(is.call(m$offset))
                  {
                    offsetFactors <- attr(Terms, "factors")
                    if(sum(offsetFactors[,offsetName]) > 1)
                      stop("Invalid argument for \"offset\": parameter of interest cannot be an interaction")
                    else
                    .offset <- eval(m$offset, envir=object$data)
                  }
                  else
                  {  
                    if(!match(offsetName, tl <- attr(Terms, "term.labels"), 
                              nomatch=FALSE))
                      .offsetName <- tl[as.logical(pmatch(tl, 
                                          offsetName, nomatch=FALSE))]
                    .offset <- do.call("eval", 
                                    list(expr=parse(text=.offsetName), 
                                         envir=object$data))
                    if(is.factor(.offset))
                    { 
                      offsetContrasts <- contrasts(.offset)
                      if(any(is.character(offsetContrasts)))
                        stop("Invalid argument for \"offset\": parameter of interest is not scalar")
                      if(any(dim(offsetContrasts) != c(2,1)))
                        stop("Invalid argument for \"offset\": parameter of interest is not scalar")
                      off.tmp <- numeric(length=length(.offset))
                      offsetLevels <- dimnames(offsetContrasts)[[1]]
                      off.tmp[.offset==offsetLevels[1]] <- 
                        offsetContrasts[1,]
                      off.tmp[.offset==offsetLevels[2]] <- 
                        offsetContrasts[2,]
                     if(any(is.na(.offset)))
                        off.tmp[is.na(.offset)] <- NA
                      .offset <- off.tmp
                    }
                  }
                })
#    assign(".offset", .offset, envir=sys.frame())                 
#    assign(".offset", .offset, pos=1)		## 20.05.13		
#    on.exit(remove(".offset", pos=1), add=TRUE)
  }
  if(!missing(pts) && (pts < 0))
    stop("Invalid argument: negative values not allowed for \"pts\"")
  else if(pts < 10)
         stop("Invalid argument: \"pts\" too small (< 10)")
  if(!missing(n) && (n < 0))
    stop("Invalid argument: negative values not allowed for \"n\"")
  else if(n < 50)
         stop("Invalid argument: \"n\" too small (< 50)")
  if(!missing(tms) && (tms < 0))
    stop("Invalid argument: negative values not allowed for \"tms\"")
  else if(tms > 1)
    warning("\"tms\" may be too large (> 1)")
  summary.obj <- summary(object)
  rsmMLE <- switch(offsetName, 
                   "scale" = log(summary.obj$dispersion[,"Estimate"]),
                   "1"     = summary.obj$coef["(Intercept)",
                                              "Estimate"],
                             summary.obj$coef[offsetName,"Estimate"])
  rsmSE <- switch(offsetName, 
                  "scale" = summary.obj$dispersion[,"Std. Error"]/
                            exp(rsmMLE),
                  "1"     = summary.obj$coef["(Intercept)",
                                             "Std. Error"],
                            summary.obj$coef[offsetName,"Std. Error"])
  wzero <- (object$weights == 0)
  anc <- object$residuals[!wzero]
  nobs <- length(anc)
  max.lik <- object$loglik
  if(offsetName != "scale")
  {
    disp <- object$dispersion
    Xmat <- model.matrix(object)[!wzero,,drop=FALSE]
    which <- match(ifelse(offsetName=="1", "(Intercept)", offsetName), 
                   dimnames(Xmat)[[2]])
    Xoffset <- Xmat[,which]
    if(!is.scalar)
    {
      if(!is.fixed)
        Xmat <- Xmat[,-which,drop=FALSE]
      tmp1 <- cbind(Xmat, anc, Xoffset)
      i.obs <- det(object$R[-which,-which,drop=FALSE]/
                   object$dispersion^2)
    }
  }
  if(!is.scalar)
  {
    dim.R <- dim(object$R)[1] - 1
    if(offsetName == "scale")
      i.obs <- det(object$R[-(dim.R+1),-(dim.R+1),drop=FALSE]/
                   object$dispersion^2)
  }
  if(is.null(from)) 
    from <- rsmMLE - 3.5*rsmSE
  if(is.null(to)) 
    to <- rsmMLE + 3.5*rsmSE
  if(from > to)
    stop("Invalid range: \"from\" < \"to\"")
  pts <- pts + (pts%%2)
  offsetCoef <- seq(from, to, length=pts)
  lengthOC <- length(offsetCoef)
  l.p <- disp.0 <- vector("numeric", lengthOC)
  if(!is.scalar) 
    part.i.obs <- vector("numeric", lengthOC)
  i.mixed <- vector("numeric", lengthOC)
  if( (offsetName != "scale") && !is.fixed && !is.scalar ) 
    part.i.mixed.1 <- vector("numeric", lengthOC)
  if((object$family[[1]] == "Huber") && object$fixed)
    object$call["dispersion"] <- object$dispersion
#  if( (length(attr(Terms, "term.labels")) == 0) &&
#      !is.null(attr(Terms, "offset")) )
#  { 
#    ..offset <- modOff <- model.offset(model.frame(object))
##    assign("..offset", ..offset, envir=sys.frame())                 
##    assign("..offset", ..offset, pos=1)		## 20.05.13		
##    on.exit(remove("..offset", pos=1), add=TRUE)
#  }
  for(i in seq(along=offsetCoef))
  {
    if(trace) 
      cat(as.name(paste("\niteration", i,
                        switch(offsetName, 
                               "scale" = ": log(scale)",
                               "1"     = ": (Intercept)", 
                                paste(": coefficient of", offsetName)),
                        "=", signif(offsetCoef[i], digits=4))))
    .object <- object$call
    .object$formula <- formula(Terms)
    switch(offsetName, 
           "scale" = { 
                       .object$dispersion <- exp(offsetCoef[i])
                       .object$data <- object$data
                     },
           "1"     = { 
                       nf <- as.formula(".~. -1")                       
                       nf <- update(.object$formula, nf)   
                       .object$formula <- nf
                       if(!is.null(.object$offset))
                         .object$offset <- object$offset + offsetCoef[i]*.offset
                       else
                         .object$offset <- offsetCoef[i]*.offset  
                       environment(.object$formula) <- 
                         environment(object$formula)             
                     },
                    { 
					   nf <- as.formula(paste(".~.-", .offsetName))			
                       nf <- update(.object$formula, nf)   
                       .object$formula <- nf
                       if(!is.null(.object$offset))
                         .object$offset <- object$offset + offsetCoef[i]*.offset
                       else
                         .object$offset <- offsetCoef[i]*.offset  
                       environment(.object$formula) <- 
                         environment(object$formula)       
                    })
#    switch(offsetName, 
#           "scale" = { 
#                       .object$dispersion <- exp(offsetCoef[i])
#                       .object$data <- object$data
#                     },
#           "1"     = { 
#                       if( (length(attr(Terms, "term.labels")) == 0) &&
#                           !is.null(attr(Terms, "offset")) )
#                       {
#                         ..offset <- modOff + offsetCoef[i]*.offset
#                         assign("..offset", ..offset, envir=sys.frame())
#                         assign("..offset", ..offset, pos=1)		## 20.05.13
#                         nf <- as.formula(
#                                 paste(deparse(formula(object)[[2]], 
#                                               width.cutoff=500), 
#                                       " ~ -1 + offset(..offset)",
#                                       collapse=""))
#                         .object$formula <- nf
#                         environment(.object$formula) <- 
#                                  environment(object$formula)
#                       }
#                       else
#                       {
#                         nf <- as.formula(
#                                 paste(deparse(formula(object),   
#                                               width.cutoff=500), 
#                                       "-1 + offset(", offsetCoef[i], 
#                                       "*.offset)", collapse=""))
#                         .object$formula <- nf
#                         environment(.object$formula) <- 
#                                  environment(object$formula)
#                       }
#                     },
#           { 
#             nf <- as.formula(paste(".~.-", .offsetName))
#             nf <- update(.object$formula, nf)
#             nf <- as.formula(paste(deparse(nf, width.cutoff=500), 
#                                    "+ offset(", offsetCoef[i], 
#                                    "* .offset)", collapse=""))
#             .object$formula <- nf
#             environment(.object$formula) <- 
#                                  environment(object$formula)
#           })
   .object <- eval(.object)
    disp.0[i] <- ifelse(offsetName == "scale", exp(offsetCoef[i]),
                        .object$dispersion[1])
    l.p[i] <- .object$loglik
    if(!is.scalar)
      part.i.obs[i] <- det(.object$R/.object$dispersion^2)  
    if(offsetName == "scale")
    {
      i.mixed[i] <- exp(rsmMLE)/disp.0[i]*sum(.object$q1[!wzero]*anc)
    }
    else if(is.fixed)
           i.mixed[i] <- sum(.object$q1[!wzero]*Xoffset)
         else 
         {
           .q2 <- .object$q2
           mat.tmp <- matrix(0, nrow=dim.R+1, ncol=dim.R+1)
           tmp2 <- cbind(Xmat, .object$residuals[!wzero])
           mat.tmp[,-1] <- t(tmp1) %*% diag(c(.q2[!wzero])) %*% tmp2
           mat.tmp[dim.R,dim.R+1] <- mat.tmp[dim.R,dim.R+1] +
                                     sum(.object$q1[!wzero]*anc)
           mat.tmp[dim.R+1,dim.R+1] <- mat.tmp[dim.R+1,dim.R+1] +
                                       sum(.object$q1[!wzero]*Xoffset)
           mat.tmp[dim.R,1] <- disp.0[i] * 
                                 (sum(.object$q1[!wzero]*anc) -
                                               nobs*disp.0[i]/disp)
           mat.tmp[dim.R+1,1] <- 
                         disp.0[i] * sum(.object$q1[!wzero]*Xoffset)
           i.mixed[i] <- det(mat.tmp/disp.0[i]^2)
           part.i.mixed.1[i] <- det(mat.tmp[1:dim.R,-1,drop=FALSE]/
                                    disp.0[i]^2)
         }
  }
  if(trace) cat("\n")
  l.mp <- l.p
  if(!is.scalar)
  {
    if(offsetName == "scale" || is.fixed)
      l.mp <- l.p - 1/2*log(part.i.obs)
    else 
      l.mp <- l.p + 1/2*log(part.i.obs) - log(part.i.mixed.1)
  }
  lp <- spline(offsetCoef[is.finite(l.p)], l.p[is.finite(l.p)], n)
  lmp <- spline(offsetCoef[is.finite(l.mp)], l.mp[is.finite(l.mp)], n)
  marg.obj <- list(workspace = list( psi=offsetCoef, l.p=lp, 
                                     l.mp=lmp ))
  marg.obj$workspace$l.p$y <- marg.obj$workspace$l.p$y -
                                max(marg.obj$workspace$l.p$y)
  marg.obj$workspace$l.mp$y <- marg.obj$workspace$l.mp$y -
                                 max(marg.obj$workspace$l.mp$y)
  if(!is.scalar)
  {
    s.mp <- predict(smooth.spline(lmp, all.knots=FALSE), lmp$x, 1)
    MLEmp <- predict(smooth.spline(s.mp$y, s.mp$x,
                                   all.knots=FALSE), 0, 0)$y
    SEmp <- sqrt(-predict(smooth.spline(s.mp$y, s.mp$x, 
                                        all.knots=FALSE), 
                          0, 1)$y)
  }
  else
  {
    MLEmp <- rsmMLE
    SEmp <- rsmSE
  }
  Diff <- rsmMLE - offsetCoef
  r.e <- Diff/rsmSE
  r.e.mp <- (MLEmp - offsetCoef)/SEmp
  r.p <- sign(Diff) * sqrt(2*(max.lik - l.p))
  if( offsetName == "scale" )
  {
    SEmp <- SEmp*exp(MLEmp)
    MLEmp <- exp(MLEmp)
  }
  MLEp <- if(offsetName == "scale") exp(rsmMLE) else rsmMLE
  SEp <- if(offsetName == "scale") rsmSE*MLEp else rsmSE
  Coef <- matrix(c(MLEp, MLEmp, SEp, SEmp), ncol=2)
  dimnames(Coef) <- list(c("uncond. ", "cond."),
                         c(" Estimate ", " Std. Error "))
  condition <- abs(offsetCoef - rsmMLE) > tms*rsmSE
  poly.ord <- round(pts/2)
  aux.mat1 <- matrix( rep(r.p,poly.ord), ncol=poly.ord ) ^
                matrix( rep(1:poly.ord,pts), ncol=poly.ord, 
                       byrow=TRUE)
  aux.mat0 <- matrix( rep(r.p,poly.ord), ncol=poly.ord ) ^
                matrix( rep(0:(poly.ord-1),pts), ncol=poly.ord, 
                       byrow=TRUE)
  if(offsetName == "scale")
  {
    rhoINF  <- SEp*(i.mixed - nobs)/MLEp
    rhoNPmp <- if(is.scalar)  rep(1, lengthOC)
               else  sqrt(part.i.obs/i.obs)
  }
  else
  {
    if(is.fixed)
    {
      rhoINF  <- SEp*i.mixed/disp.0
      rhoNPmp <- if(is.scalar)  rep(1, lengthOC)
                 else  sqrt(part.i.obs/i.obs)
    }
    else 
    {
      rhoINF  <- SEp*i.mixed/part.i.mixed.1
      rhoNPmp <- if(is.scalar)  rep(1, lengthOC)
                 else part.i.mixed.1/sqrt(i.obs*part.i.obs)
    }
  }
  rhoLR <- rhoINF*rhoNPmp
  rhoLR <- ifelse(rhoLR*r.p < 0, -rhoLR, rhoLR)
  q.mp <- spline(x = offsetCoef[is.finite(rhoLR)], 
                 y = rhoLR[is.finite(rhoLR)], n)
  aux.mod <- lm.fit(aux.mat1, rhoINF)
  log.rhoINF <- ifelse(condition, log(abs(rhoINF/r.p)),
                                  log(abs(aux.mat0%*%aux.mod$coef)))
  log.rho <- log(abs(rhoNPmp)) + log.rhoINF
##
##-- INF 
##
  INFmp <- INFmp.rp <- log.rhoINF
##
##-- NP
##
  if(!is.scalar)
    NPmp <- NPmp.rp <- log(abs(rhoNPmp)) 
  else 
    NPmp <- NPmp.rp <- rep(0, length=lengthOC)
##
## step 2: log(\rho_INF/r.p)/r.p + log(\rho_NP)/r.p  
##
##-- r*
##
  aux.mod <- lm.fit(aux.mat1, log.rho)
  rho <- aux.mat0%*%aux.mod$coef
  rho <- ifelse(condition, log.rho/r.p, rho)
  r.mp <- r.p + rho
##
##-- INF
##
  INFmp.tmp <- INFmp
  aux.mod <- lm.fit(aux.mat1, INFmp.tmp)
  INFmp.tmp <- aux.mat0%*%aux.mod$coef
  INFmp <- ifelse(condition, INFmp/r.p, INFmp.tmp)
##
##-- NP
##
  if(!is.scalar)
  {
    NPmp.tmp <- NPmp
    aux.mod <- lm.fit(aux.mat1, NPmp.tmp)
    NPmp.tmp <- aux.mat0%*%aux.mod$coef
    NPmp <- ifelse(condition, NPmp/r.p, NPmp.tmp)
  }
##
##-- L.-R.
##
  rhoLR.tmp <- (1/rhoLR - 1/r.p)[condition]
  aux.mod <- lm.fit(aux.mat0[condition,], rhoLR.tmp)
  rhoLR.tmp <- aux.mat0%*%aux.mod$coef
  rhoLR <- ifelse(condition, 1/rhoLR - 1/r.p, rhoLR.tmp)
  LR <- pnorm(r.p) - dnorm(r.p)*rhoLR
##
##--> stop
##
  marg.obj$workspace <-
          c(marg.obj$workspace,
            list( r.e        = 
                          (re <- spline(offsetCoef[is.finite(r.e)], 
                                        r.e[is.finite(r.e)], n)),
                  r.e.mp     = if(!is.scalar)
                                 spline(offsetCoef[is.finite(r.e.mp)],
                                        r.e.mp[is.finite(r.e.mp)], n)
                               else re,
                  r.p        = spline(offsetCoef[is.finite(r.p)], 
                                      r.p[is.finite(r.p)], n),
                  r.mp       = spline(offsetCoef[is.finite(r.mp)], 
                                      r.mp[is.finite(r.mp)], n) ))
  lr.mp <- spline(offsetCoef[is.finite(LR)], LR[is.finite(LR)], n)
  lr.mp$y[lr.mp$y < 0] <- 0
  lr.mp$y[lr.mp$y > 1] <- 1
  marg.obj$workspace <-
          c(marg.obj$workspace,
            list( lr.mp      = lr.mp, 
                  q          = q.mp, 
                  inf.mp     = spline(offsetCoef[is.finite(INFmp)], 
                                      INFmp[is.finite(INFmp)], n),
                  inf.mp.rp  = spline(r.p[is.finite(INFmp.rp)], 
                                      INFmp.rp[is.finite(INFmp.rp)], 
                                      n),
                  np.mp      = if(!is.scalar)
                                 spline(offsetCoef[is.finite(NPmp)], 
                                        NPmp[is.finite(NPmp)], n)
                               else NULL,
                  np.mp.rp   = if(!is.scalar) 
                                 spline(r.p[is.finite(NPmp.rp)], 
                                        NPmp.rp[is.finite(NPmp.rp)], 
                                        n)
                               else NULL
                ))
  marg.obj <- c(marg.obj, 
                list( coefficients = Coef,
                      call         = m,
                      formula      = object$formula,
                      family       = as.name(object$dist),
                      offset       = as.name(offsetName),
                      diagnostics  = c(INF = 
                          max(abs(marg.obj$workspace$inf.mp$y)),
                                       NP  = if(!is.scalar)
                          max(abs(marg.obj$workspace$np.mp$y))
                                               else NULL),
                      n.approx     = pts,
                      omitted.val  = c(rsmMLE - tms*rsmSE,
                                       rsmMLE + tms*rsmSE),
                      is.scalar    = is.scalar ))
  attr(marg.obj, "class") <- c("marg", "cond")
  marg.obj
}


print.marg <- function(x, digits = max(3, getOption("digits")-3), ...)
{
  is.scalar <- x$is.scalar
  if(is.scalar)
    cat("Original model contains no nuisance parameter\n\n")
  if(!is.null(cl <- x$call)) 
  {
    cat("Call:\n")
    dput(cl)
  }
  cat("\nFormula:  ")
  dput(x$formula)
  cat("Family:  ")
  print(x$family)
  cat("Offset:  ")
  if( x$offset != "1" )  
    print(as.name(x$offset))
  else cat("Intercept\n")
  cat("\n")
  if( is.scalar )
    print(x$coef[1,], digits=digits)
  else print(x$coef, digits=digits)
  xd <- x$diagnostics
  cat("\nDiagnostics: \n")
  if( is.scalar )
    print(xd[1], digits=digits)
  else print(xd, digits=digits)
  cat("\n Approximation based on", x$n.approx, "points\n")
}

plot.marg <- function(x = stop("nothing to plot"), from = x.axis[1],
                      to = x.axis[n], which = NULL, alpha = 0.05,
                      add.leg = TRUE, loc.leg = FALSE, 
                      add.labs = TRUE, cex = 0.7, cex.lab = 1, 
                      cex.axis = 1, cex.main = 1, lwd1 = 1, lwd2 = 2,
                      lty1 = "solid", lty2 = "dashed", 
                      col1 = "black", col2 = "blue", 
                      tck = 0.02, las = 1, adj = 0.5, 
                      lab = c(15, 15, 5), ...)             
{
  offset <- x$offset
  if( offset == "1" )
    offset <- "Intercept"
  is.scalar <- x$is.scalar
  if(is.scalar)
    cat("\nOriginal model contains no nuisance parameter\n")
  if( x$family == "Huber" )
     cat(paste("\nNOTE: The following plots are only meaningful",
               " for Huber's least favourable distribution\n",
               "if the errors are _truly_",
               " Huber-type distributed.\n", sep=""))
  choices <- c("All",
               ifelse(!is.scalar, 
 	         "Profile and modified profile log likelihoods",
 	         "Profile log likelihood"),
               ifelse(!is.scalar, 
  	         "Profile and modified profile likelihood ratios",
                 "Profile likelihood ratio"),
	       "Profile and modified likelihood roots",
               "Lugannani-Rice approximation",
               "Confidence intervals",
               ifelse(!is.scalar,
                 "Diagnostics based on INF/NP decomposition\n",
                 "Diagnostics based on INF decomposition\n"))
  tmenu <- paste("", choices)
  if( is.null(which) )
    pick <- menu(tmenu, 
                 title = "\n Make a plot selection (or 0 to exit)\n")
  else if( !match(which, 2:7, nomatch=FALSE) )
               stop("choice not valid")
       else pick <- which
  if(pick == 0)
    stop(" no graph required ! ")
#  attach(x$workspace, warn.conflicts = FALSE)
#  on.exit(invisible(detach()))
  r.p <- x$workspace$r.p
  r.mp <- x$workspace$r.mp
  l.p <- x$workspace$l.p
  l.mp <- x$workspace$l.mp
  r.e <- x$workspace$r.e
  r.e.mp <- x$workspace$r.e.mp
  lr.mp <- x$workspace$lr.mp
  inf.mp <- x$workspace$inf.mp
  inf.mp.rp <- x$workspace$inf.mp.rp
  np.mp <- x$workspace$np.mp
  np.mp.rp <- x$workspace$np.mp.rp
  is.scale <- (paste(x$offset) == "scale")
  coeff <- x$coefficients
  if(is.scale)
  {
    coeff[,2] <- coeff[,2]/coeff[,1]
    coeff[,1] <- log(coeff[,1])
  }
  x.alpha <- qnorm(1 - alpha/2)
  CI.rp <- predict(smooth.spline(r.p$y, r.p$x), 
                   c(x.alpha, -x.alpha), 0)$y
  CI.rmp <- predict(smooth.spline(r.mp$y, r.mp$x), 
                    c(x.alpha, -x.alpha), 0)$y
  CI.mle <- (coeff[1,1] - coeff[1,2]*c(x.alpha, -x.alpha))
  CI.cmle <- (coeff[2,1] - coeff[2,2]*c(x.alpha, -x.alpha))
  if(is.scale)
  {
    CI.rp <- exp(CI.rp)
    CI.rmp <- exp(CI.rmp)
    CI.mle <- exp(CI.mle)
    CI.cmle <- exp(CI.cmle)
  }
  from.ic <- min(c(CI.rp[1], CI.rmp[1], CI.mle[1], CI.cmle[1]))
  to.ic <- max(c(CI.rp[2], CI.rmp[2], CI.mle[2], CI.cmle[2]))
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par), add=TRUE)
  repeat
  {
    invisible(par(tck=tck, las=las, adj=adj, lab=lab, ...))
    switch(pick,
		"1" = {
			invisible(split.screen(figs = c(1,1)))
##               ---------------
           screen(1)
           n <- length(l.p$x)
           x.axis <- l.p$x  
           condition <- (x.axis >= from) & (x.axis <= to) &
                          ((l.p$y > -3) | (l.mp$y > -3)) 
           plot(x.axis[condition], l.p$y[condition], 
                type = "n", lty = lty1, lwd=lwd1, 
                ylim = c(max(-3, 
                         min(l.p$y[condition], l.mp$y[condition])), 0),
                xlab = "", ylab= "", cex.lab = cex.lab, 
                cex.axis=cex.axis, cex.main=cex.main, ...)
           lines(x.axis[condition], l.p$y[condition], lty=lty1, 
                 col=col1, lwd=lwd1, ...)
           if( !is.scalar )
             lines(x.axis[condition], l.mp$y[condition], lty=lty1, 
                   col=col2, lwd=lwd2, ...)
           if(add.labs)
              title(xlab = if(is.scale) "log(scale)" else
                             paste("coefficient of", offset),
                    ylab = "log likelihood",
                    main = ifelse(!is.scalar, 
                                  "Profile and modified profile log LIKs",
                                  "Profile log likelihood"),
                    ...)
	   if(add.leg == TRUE)
           {
             if(loc.leg)
             {
               cat("Choose legend position\n")
               loci <- locator(1)
               x.lim <- loci$x
               y.lim <- loci$y
             }
             else
             {
               x.lim <- (x.axis[l.mp$y == max(l.mp$y)] +
                           x.axis[condition][1])/2
               y.lim <- max(-3,
                            min(l.p$y[condition], l.mp$y[condition]))/2
             }
             if( !is.scalar )
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("profile log likelihood",
                        "modified profile log likelihood"),
                      cex=cex, lty=c(lty1,lty1), lwd=c(lwd1,lwd2), 
                      col=c(col1,col2), bty="n", ...)
             else
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      "profile log likelihood",
                      cex=cex, lty=lty1, lwd=lwd1, 
                      col=col1, bty="n", ...)

           }
##               ---------------
                        invisible(close.screen(all.screens = TRUE))
                        par(ask = TRUE)
                        invisible(split.screen(figs = c(2,2)))
##               ---------------
           screen(1)
           conf.limit <-  qchisq(1-alpha, 1)
           condition <- (x.axis >= from) & (x.axis <= to) &
                        ( ((r.e$y)^2 < (conf.limit + 1)) |
                          ((r.e.mp$y)^2 < (conf.limit + 1)) |
                          (-2 * l.p$y < (conf.limit + 1)) |
                          (-2 * l.mp$y < (conf.limit + 1)) )
           if( !is.scalar )
           {
             matplot(x.axis[condition],
                     cbind( ((r.e$y)^2)[condition], 
                            ((r.e.mp$y)^2)[condition],
                            (-2 * l.p$y)[condition] ),
                     xlab = "", ylab= "", type = "l", 
                     lty = c(lty2,lty2,lty1), lwd=c(lwd1,lwd2,lwd1), 
                     col=c(col1,col2,col1), cex.lab=cex.lab,
                     cex.main=cex.main, cex.axis=cex.axis, ...)
  	     lines(x.axis[condition], (-2 * l.mp$y)[condition], 
                   lty=lty1, lwd=lwd2, col=col2, ...)
           }
           else
             matplot(x.axis[condition],
                     cbind( ((r.e$y)^2)[condition], 
                            (-2 * l.p$y)[condition] ),
                     xlab = "", ylab= "", type = "l", 
                     lty = c(lty2,lty1), 
                     lwd=c(lwd1,lwd1), 
                     col=c(col1,col1), cex.lab=cex.lab,
                     cex.main=cex.main, cex.axis=cex.axis, ...)
           abline(h=conf.limit, lty="dotted", ...)
           if(add.labs)
             title(xlab = if(is.scale) "log(scale)" else
                            paste("coefficient of", offset),
                   ylab = "likelihood ratio",
                   main = ifelse(!is.scalar,
                                 "Profile and modified profile LRs",
                                 "Profile likelihood ratio"),
                   ...)
	   if(add.leg)
           {
             if(loc.leg)
             {
               cat("Choose legend position\n")
               loci <- locator(1)
               x.lim <- loci$x
               y.lim <- loci$y
             }
             else
             {
               x.lim <- ( x.axis[l.p$y == max(l.p$y)] +
                          x.axis[condition][1] )/2
               tt <- max(c(max(((r.e$y)^2)[condition]),
                           max(((r.e.mp$y)^2)[condition]),
                           max((-2 * l.p$y)[condition]),
                           max((-2 * l.mp$y)[condition])))
	       y.lim <- if(tt > conf.limit) (tt+conf.limit)/1.9
                          else tt/2
             }
             if( !is.scalar )
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("profile likelihood ratio", "modified profile LR",
                        "Wald pivot",
                        "Wald pivot (cond.)"),
                      lty=c(lty1,lty1,lty2,lty2), 
                      lwd=c(lwd1,lwd2,lwd1,lwd2),
                      col=c(col1,col2,col1,col2), bty="n", 
                      cex=cex, ...)
             else 
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("profile likelihood ratio", "Wald pivot"),
                      lty=c(lty1,lty2), lwd=c(lwd1,lwd1),
                      col=c(col1,col1), bty="n", cex=cex, ...)
             }
##               ---------------
           conf.limit <- qnorm(1 - alpha/2)
           screen(2)
           condition <- (x.axis >= from) & (x.axis <= to) &
                          ( ((r.e$y < 4) & (r.e$y > -4)) |
                            ((r.e.mp$y < 4) & (r.e.mp$y > -4)) |
                            ((r.p$y < 4) & (r.p$y > -4)) |
                            ((r.mp$y < 4) & (r.mp$y > -4)) )
           plot(0, 0, type="n", xlim=range(x.axis[condition]), ylim=c(-4, 4),
                xlab="", ylab="", cex.lab=cex.lab, cex.axis=cex.axis, 
                cex.main=cex.main, ...)
#           condition <- (x.axis >= from) & (x.axis <= to) &
#                          (r.e$y < 4) & (r.e$y > -4)
	   lines(x.axis[condition], r.e$y[condition], lty=lty2, 
                 lwd=lwd1, col=col1, ...)
           if( !is.scalar )
           {
#             condition <- (x.axis >= from) & (x.axis <= to) &
#                            (r.e.mp$y < 4) & (r.e.mp$y > -4)
	     lines(x.axis[condition], r.e.mp$y[condition], lty=lty2, 
                   lwd=lwd2, col=col2, ...)
           }
#           condition <- (x.axis >= from) & (x.axis <= to) &
#                          (r.p$y < 4) & (r.p$y > -4)
	   lines(x.axis[condition], r.p$y[condition], lty=lty1, 
                 lwd=lwd1, col=col1, ...)
#           condition <- (x.axis >= from) & (x.axis <= to) &
#                          (r.mp$y < 4) & (r.mp$y > -4)
	   lines(x.axis[condition], r.mp$y[condition], lty=lty1,
                 lwd=lwd2, col=col2, ... )
	   abline(h=0, lty="dotted", ...)
	   abline(h=conf.limit, lty="dotted", ...)
	   abline(h=-conf.limit, lty="dotted", ...)
           if(add.labs)
             title(xlab = if(is.scale) "log(scale)" else
                            paste("coefficient of", offset),
                   ylab = "likelihood root",
                   main = "Profile and modified likelihood roots", 
                   ...)
	   if(add.leg == TRUE)
           {
             if(loc.leg)
             {
               cat("Choose legend position\n")
               loci <- locator(1)
               x.lim <- loci$x
               y.lim <- loci$y
             }
             else
             {
               x.lim <- x.axis[condition][2]
               y.lim <- -2
             }
             if( !is.scalar ) 
  	       legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("likelihood root", 
                        "modified likelihood root",
                        "Wald pivot", 
                        "Wald pivot (cond.)"),
                       lty=c(lty1,lty1,lty2,lty2), 
                       lwd=c(lwd1,lwd2,lwd1,lwd2), 
                       col=c(col1,col2,col1,col2), bty="n", 
                       cex=cex, ...)
             else
  	       legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("likelihood root", 
                        "modified likelihood root",
                        "Wald pivot"),
                       lty=c(lty1,lty1,lty2), 
                       lwd=c(lwd1,lwd2,lwd1), 
                       col=c(col1,col2,col1), bty="n", 
                       cex=cex, ...) 
           }
##               ---------------
           screen(3)
           n <- length(r.mp$x)
           x.axis <- r.mp$x  
           condition <- (x.axis >= from) & (x.axis <= to)
           plot(x.axis[condition], lr.mp$y[condition], ylim=c(0,1),
                xlab="", ylab="", type="l", lty=lty1, lwd=lwd2,
                col=col2, cex.lab=cex.lab, cex.axis=cex.axis, 
                cex.main=cex.main, ...)
           if(add.labs)
             title(xlab = if(is.scale) "log(scale)" else
                            paste("coefficient of", offset),
                   ylab = "tail area approximation",
                   main = "Lugannani-Rice approximation",
                   cex = cex, ...)
##               ---------------
           screen(4)
           invisible(par(lab = c(30,5,5)))
           plot(exp(coeff[1,1]), 3.5, type="n",
                xlim=c(from.ic, to.ic), ylim=c(0,ifelse(!is.scalar, 5, 4)), 
                xlab="", ylab="", yaxt="n", cex.lab=cex.lab, 
                cex.main=cex.main, cex.axis=cex.axis, ...)
           segments(CI.rp[1], 2, CI.rp[2], 2, lty=lty1, lwd=lwd1, 
                    col=col1, ...)
           segments(CI.rmp[1], 1, CI.rmp[2], 1, lty=lty1, lwd=lwd1, 
                    col=col2, ...)
           if( !is.scalar )
           {
             segments(CI.cmle[1], 3, CI.cmle[2], 3, lty=lty1, lwd=lwd1,
                      col=col2, ...)
             segments(CI.mle[1], 4, CI.mle[2], 4, lty=lty1, lwd=lwd1, 
                      col=col1, ...)
           }
           else
             segments(CI.mle[1], 3, CI.mle[2], 3, lty=lty1, lwd=lwd1, 
                      col=col1, ...)
	   if(add.labs)
             title(xlab = if(is.scale) "scale" else
                            paste("coefficient of", offset),
                   ylab = "confidence interval",
                   main = paste(100*(1-alpha), 
                                "% Confidence Intervals "),
                   ...)
           if(add.leg == TRUE)
           {
             text(if(is.scale) exp(coeff[2,1]) else coeff[2,1], 2.2,
                  "likelihood root", cex=cex, ...)
             text(if(is.scale) exp(coeff[2,1]) else coeff[2,1], 1.2,
                  "modified likelihood root",
                  cex=cex, ...)
             if( !is.scalar )
             {
               text(if(is.scale) exp(coeff[2,1]) else coeff[2,1], 4.2,
                    "Wald pivot",
                    cex=cex, ...)
               text(if(is.scale) exp(coeff[2,1]) else coeff[2,1], 3.2,
                    "Wald pivot (cond.)",
                    cex=cex, ...)
             }
             else
               text(if(is.scale) exp(coeff[2,1]) else coeff[2,1], 3.2,
                    "Wald pivot",
                    cex=cex, ...)
           }
##               ---------------
			invisible(close.screen(all.screens = TRUE))
                        par(ask = TRUE)
                        par(pty="s")
                        if(is.scalar)  split.screen(figs=c(1,2))
                          else  split.screen(figs=c(2,2))
##               ---------------
           screen(1)
           plot(inf.mp, type="n", xlab="", ylab="", cex.lab=cex.lab,
                cex.axis=cex.axis, cex.main=cex.main, ...)
           lines(inf.mp, lty=lty1, lwd=lwd1, col=col1, ...)
           abline(h=0.2, lty=3, ...)
           abline(h=-0.2, lty=3, ...)
           if(add.labs)
             title(xlab=if(is.scale) "log(scale)" else
                          paste("coefficient of", offset), ylab="INF",
                   main=" INF correction term ", ...)
           screen(2)
           plot(inf.mp.rp, type="l", xlab="", ylab="", lty=lty1, 
                lwd=lwd1, cex.lab=cex.lab, cex.axis=cex.axis, 
                cex.main=cex.main, col=col1, ...)
           if(add.labs)
             title(xlab="likelihood root", ylab="INF",
                   main="INF correction term", ...)
           if(!is.scalar)
           {
             screen(3)
             plot(np.mp, type="n", xlab="", ylab="", cex.lab=cex.lab,
                  cex.axis=cex.axis, cex.main=cex.main, ...)
             lines(np.mp, lty=lty1, lwd=lwd1, col=col1, ...)
             abline(h=0.2, lty=3, ...)
             abline(h=-0.2, lty=3, ...)
             if(add.labs)
                title(xlab=if(is.scale) "log(scale)" else
                             paste("coefficient of", offset), ylab="NP",
                      main=" NP correction term ", ...)
             screen(4)
             plot(np.mp.rp, type="l", xlab="", ylab="", lty=lty1,
                  lwd=lwd1, cex.lab=cex.lab, cex.axis=cex.axis,
                  cex.main=cex.main, col=col1, ...)
             if(add.labs)
               title(xlab="likelihood root", ylab="NP",
                     main="NP correction term", ...)
           }
##               ---------------
                        close.screen(all.screens=TRUE)
                        par(pty="m")
                        par(ask = FALSE)
		}
		,
		"2" = {
                        invisible(par(mfrow=c(1,1)))
##               ---------------
           n <- length(l.p$x)
           x.axis <- l.p$x  
           condition <- (x.axis >= from) & (x.axis <= to) &
                          ((l.p$y > -3) | (l.mp$y > -3)) 
           plot(x.axis[condition], l.p$y[condition], 
                type = "n", lty = lty1, lwd=lwd1, 
                ylim = c(max(-3, 
                         min(l.p$y[condition], l.mp$y[condition])), 0),
                xlab = "", ylab= "", cex.lab = cex.lab, 
                cex.axis=cex.axis, cex.main=cex.main, ...)
           lines(x.axis[condition], l.p$y[condition], lty=lty1, 
                 col=col1, lwd=lwd1, ...)
           if( !is.scalar )
             lines(x.axis[condition], l.mp$y[condition], lty=lty1, 
                   col=col2, lwd=lwd2, ...)
           if(add.labs)
              title(xlab = if(is.scale) "log(scale)" else
                             paste("coefficient of", offset),
                    ylab = "log likelihood",
                    main = ifelse(!is.scalar, 
                                  "Profile and modified profile log LIKs",
                                  "Profile log likelihood"),
                    ...)
	   if(add.leg == TRUE)
           {
             if(loc.leg)
             {
               cat("Choose legend position\n")
               loci <- locator(1)
               x.lim <- loci$x
               y.lim <- loci$y
             }
             else
             {
               x.lim <- (x.axis[l.mp$y == max(l.mp$y)] +
                           x.axis[condition][1])/2
               y.lim <- max(-3,
                            min(l.p$y[condition], l.mp$y[condition]))/2
             }
             if( !is.scalar )
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("profile log likelihood",
                        "modified profile log likelihood"),
                      cex=cex, lty=c(lty1,lty1), lwd=c(lwd1,lwd2), 
                      col=c(col1,col2), bty="n", ...)
             else
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      "profile log likelihood",
                      cex=cex, lty=lty1, lwd=lwd1, 
                      col=col1, bty="n", ...)
           }  
##               ---------------
                }
		,
		"3" = {
                        invisible(par(mfrow=c(1,1)))
##               ---------------
           n <- length(l.p$x)
           x.axis <- l.p$x  
           conf.limit <-  qchisq(1-alpha, 1)
           condition <- (x.axis >= from) & (x.axis <= to) &
                        ( ((r.e$y)^2 < (conf.limit + 1)) |
                          ((r.e.mp$y)^2 < (conf.limit + 1)) |
                          (-2 * l.p$y < (conf.limit + 1)) |
                          (-2 * l.mp$y < (conf.limit + 1)) )
           if( !is.scalar )
           {
             matplot(x.axis[condition],
                     cbind( ((r.e$y)^2)[condition], 
                            ((r.e.mp$y)^2)[condition],
                            (-2 * l.p$y)[condition] ),
                     xlab = "", ylab= "", type = "l", 
                     lty = c(lty2,lty2,lty1), lwd=c(lwd1,lwd2,lwd1), 
                     col=c(col1,col2,col1), cex.lab=cex.lab,
                     cex.main=cex.main, cex.axis=cex.axis, ...)
  	     lines(x.axis[condition], (-2 * l.mp$y)[condition], 
                   lty=lty1, lwd=lwd2, col=col2, ...)
           }
           else
             matplot(x.axis[condition],
                     cbind( ((r.e$y)^2)[condition], 
                            (-2 * l.p$y)[condition] ),
                     xlab = "", ylab= "", type = "l", 
                     lty = c(lty2,lty1), 
                     lwd=c(lwd1,lwd1), 
                     col=c(col1,col1), cex.lab=cex.lab,
                     cex.main=cex.main, cex.axis=cex.axis, ...)
           abline(h=conf.limit, lty="dotted", ...)
           if(add.labs)
             title(xlab = if(is.scale) "log(scale)" else
                            paste("coefficient of", offset),
                   ylab = "likelihood ratio",
                   main = ifelse(!is.scalar,
                                 "Profile and modified profile LRs",
                                 "Profile likelihood ratio"),
                   ...)
	   if(add.leg)
           {
             if(loc.leg)
             {
               cat("Choose legend position\n")
               loci <- locator(1)
               x.lim <- loci$x
               y.lim <- loci$y
             }
             else
             {
               x.lim <- ( x.axis[l.p$y == max(l.p$y)] +
                          x.axis[condition][1] )/2
               tt <- max(c(max(((r.e$y)^2)[condition]),
                           max(((r.e.mp$y)^2)[condition]),
                           max((-2 * l.p$y)[condition]),
                           max((-2 * l.mp$y)[condition])))
	       y.lim <- if(tt > conf.limit) (tt+conf.limit)/1.9
                          else tt/2
             }
             if( !is.scalar )
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("profile likelihood ratio", "modified profile LR",
                        "Wald pivot",
                        "Wald pivot (cond.)"),
                      lty=c(lty1,lty1,lty2,lty2), 
                      lwd=c(lwd1,lwd2,lwd1,lwd2),
                      col=c(col1,col2,col1,col2), bty="n", 
                      cex=cex, ...)
             else 
               legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("profile likelihood ratio", "Wald pivot"),
                      lty=c(lty1,lty2), lwd=c(lwd1,lwd1),
                      col=c(col1,col1), bty="n", cex=cex, ...)
           }
##               ---------------
       	        }
		,
		"4" = {
                        invisible(par(mfrow=c(1,1)))
##               ---------------
           n <- length(r.mp$x)
           x.axis <- r.mp$x  
           conf.limit <-  qnorm(1-alpha/2)
           condition <- (x.axis >= from) & (x.axis <= to) &
                          ( ((r.e$y < 4) & (r.e$y > -4)) |
                            ((r.e.mp$y < 4) & (r.e.mp$y > -4)) |
                            ((r.p$y < 4) & (r.p$y > -4)) |
                            ((r.mp$y < 4) & (r.mp$y > -4)) )
           plot(0, 0, type="n", xlim=range(x.axis[condition]), ylim=c(-4, 4),
                xlab="", ylab="", cex.lab=cex.lab, cex.axis=cex.axis, 
                cex.main=cex.main, ...)
#           condition <- (x.axis >= from) & (x.axis <= to) &
#                          (r.e$y < 4) & (r.e$y > -4)
	   lines(x.axis[condition], r.e$y[condition], lty=lty2, 
                 lwd=lwd1, col=col1, ...)
           if( !is.scalar )
           {
#             condition <- (x.axis >= from) & (x.axis <= to) &
#                            (r.e.mp$y < 4) & (r.e.mp$y > -4)
	     lines(x.axis[condition], r.e.mp$y[condition], lty=lty2, 
                   lwd=lwd2, col=col2, ...)
           }
#           condition <- (x.axis >= from) & (x.axis <= to) &
#                          (r.p$y < 4) & (r.p$y > -4)
	   lines(x.axis[condition], r.p$y[condition], lty=lty1, 
                 lwd=lwd1, col=col1, ...)
#           condition <- (x.axis >= from) & (x.axis <= to) &
#                          (r.mp$y < 4) & (r.mp$y > -4)
	   lines(x.axis[condition], r.mp$y[condition], lty=lty1,
                 lwd=lwd2, col=col2, ... )
	   abline(h=0, lty="dotted", ...)
	   abline(h=conf.limit, lty="dotted", ...)
	   abline(h=-conf.limit, lty="dotted", ...)
           if(add.labs)
             title(xlab = if(is.scale) "log(scale)" else
                            paste("coefficient of", offset),
                   ylab = "likelihood root",
                   main = "Profile and modified likelihood roots", 
                   ...)
	   if(add.leg == TRUE)
           {
             if(loc.leg)
             {
               cat("Choose legend position\n")
               loci <- locator(1)
               x.lim <- loci$x
               y.lim <- loci$y
             }
             else
             {
               x.lim <- x.axis[condition][2]
               y.lim <- -2
             }
             if( !is.scalar ) 
  	       legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("likelihood root", 
                        "modified likelihood root",
                        "Wald pivot", 
                        "Wald pivot (cond.)"),
                       lty=c(lty1,lty1,lty2,lty2), 
                       lwd=c(lwd1,lwd2,lwd1,lwd2), 
                       col=c(col1,col2,col1,col2), bty="n", 
                       cex=cex, ...)
             else
  	       legend(c(x.lim, x.lim), c(y.lim, y.lim),
                      c("likelihood root", 
                        "modified likelihood root",
                        "Wald pivot"),
                       lty=c(lty1,lty1,lty2), 
                       lwd=c(lwd1,lwd2,lwd1), 
                       col=c(col1,col2,col1), bty="n", 
                       cex=cex, ...) 
           }
##               ---------------
		}
		,
		"5" = {
                        invisible(par(mfrow=c(1,1)))
##               ---------------
           n <- length(r.mp$x)
           x.axis <- r.mp$x  
           condition <- (x.axis >= from) & (x.axis <= to)
           plot(x.axis[condition], lr.mp$y[condition], ylim=c(0,1),
                xlab="", ylab="", type="l", lty=lty1, lwd=lwd2,
                col=col2, cex.lab=cex.lab, cex.axis=cex.axis, 
                cex.main=cex.main, ...)
           if(add.labs)
             title(xlab = if(is.scale) "log(scale)" else
                            paste("coefficient of", offset),
                   ylab = "tail area approximation",
                   main = "Lugannani-Rice approximation",
                   cex = cex, ...)
                }
##               ---------------
                ,
                "6" = {
                        invisible(par(lab = c(30, 5, 5)))
           plot(exp(coeff[1,1]), 3.5, type="n",
                xlim=c(from.ic, to.ic), ylim=c(0,ifelse(!is.scalar, 5, 4)), 
                xlab="", ylab="", yaxt="n", cex.lab=cex.lab, 
                cex.main=cex.main, cex.axis=cex.axis, ...)
           segments(CI.rp[1], 2, CI.rp[2], 2, lty=lty1, lwd=lwd1, 
                    col=col1, ...)
           segments(CI.rmp[1], 1, CI.rmp[2], 1, lty=lty1, lwd=lwd1, 
                    col=col2, ...)
           if( !is.scalar )
           {
             segments(CI.cmle[1], 3, CI.cmle[2], 3, lty=lty1, lwd=lwd1,
                      col=col2, ...)
             segments(CI.mle[1], 4, CI.mle[2], 4, lty=lty1, lwd=lwd1, 
                      col=col1, ...)
           }
           else
             segments(CI.mle[1], 3, CI.mle[2], 3, lty=lty1, lwd=lwd1, 
                      col=col1, ...)
	   if(add.labs)
             title(xlab = if(is.scale) "scale" else
                            paste("coefficient of", offset),
                   ylab = "confidence interval",
                   main = paste(100*(1-alpha), 
                                "% Confidence Intervals "),
                   ...)
           if(add.leg == TRUE)
           {
             text(if(is.scale) exp(coeff[2,1]) else coeff[2,1], 2.2,
                  "likelihood root", cex=cex, ...)
             text(if(is.scale) exp(coeff[2,1]) else coeff[2,1], 1.2,
                  "modified likelihood root",
                  cex=cex, ...)
             if( !is.scalar )
             {
               text(if(is.scale) exp(coeff[2,1]) else coeff[2,1], 4.2,
                    "Wald pivot",
                    cex=cex, ...)
               text(if(is.scale) exp(coeff[2,1]) else coeff[2,1], 3.2,
                    "Wald pivot (cond.)",
                    cex=cex, ...)
             }
             else
               text(if(is.scale) exp(coeff[2,1]) else coeff[2,1], 3.2,
                    "Wald pivot",
                    cex=cex, ...)
           }
##               ---------------
                }
                , 
                "7" = {
			invisible(close.screen(all.screens = TRUE))
                        par(pty="s")
                        if(is.scalar)  split.screen(figs=c(1,2))
                          else  split.screen(figs=c(2,2))
##               ---------------

           screen(1)
           plot(inf.mp, type="n", xlab="", ylab="", cex.lab=cex.lab,
                cex.axis=cex.axis, cex.main=cex.main, ...)
           lines(inf.mp, lty=lty1, lwd=lwd1, col=col1, ...)
           abline(h=0.2, lty=3, ...)
           abline(h=-0.2, lty=3, ...)
           if(add.labs)
             title(xlab=if(is.scale) "log(scale)" else
                          paste("coefficient of", offset), ylab="INF",
                   main=" INF correction term ", ...)
           screen(2)
           plot(inf.mp.rp, type="l", xlab="", ylab="", lty=lty1, 
                lwd=lwd1, cex.lab=cex.lab, cex.axis=cex.axis, 
                cex.main=cex.main, col=col1, ...)
           if(add.labs)
             title(xlab="likelihood root", ylab="INF",
                   main="INF correction term", ...)
           if(!is.scalar)
           {
             screen(3)
             plot(np.mp, type="n", xlab="", ylab="", cex.lab=cex.lab,
                  cex.axis=cex.axis, cex.main=cex.main, ...)
             lines(np.mp, lty=lty1, lwd=lwd1, col=col1, ...)
             abline(h=0.2, lty=3, ...)
             abline(h=-0.2, lty=3, ...)
             if(add.labs)
                title(xlab=if(is.scale) "log(scale)" else
                             paste("coefficient of", offset), ylab="NP",
                      main=" NP correction term ", ...)
             screen(4)
             plot(np.mp.rp, type="l", xlab="", ylab="", lty=lty1,
                  lwd=lwd1, cex.lab=cex.lab, cex.axis=cex.axis,
                  cex.main=cex.main, col=col1, ...)
             if(add.labs)
               title(xlab="likelihood root", ylab="NP",
                     main="NP correction term", ...)
           }
##               ---------------
                        close.screen(all.screens=TRUE)
                        par(pty="m")
                        })
    if( missing(which) )
      pick <- menu(tmenu, 
                   title="\n Make a plot selection (or 0 to exit)\n")
    if( (pick == 0) || !missing(which) )
    {
      invisible(close.screen(all.screens=TRUE))
      break
    }
  }
}

summary.marg <- function(object, alpha = 0.05, test = NULL, 
                         all = FALSE, coef = TRUE, 
                         int = ifelse((is.null(test) || all), 
                                      TRUE, FALSE),
                         digits = NULL, ...)
{
  m <- match.call()
  dim.alpha <- length(alpha)
  if( !is.null(test) )
    dim.test <- length(test)
  alpha.quant <- c(qnorm(1 - alpha/2), qnorm(1 - alpha))
#  attach(object$workspace, warn.conflicts = FALSE)
#  on.exit( detach() )
  r.p <- object$workspace$r.p
  r.mp <- object$workspace$r.mp
  lr.mp <- object$workspace$lr.mp
  cf <- object$coefficients
  is.scale <- ( paste(object$offset) == "scale" )
  is.scalar <- object$is.scalar
  if(int)
  {
    if(is.scale)
    {
              cf[,2] <- cf[,2] / cf[,1]
              cf[,1] <- log(cf[,1])
    }
    CI.mle <- cf[1,1] - cf[1,2] * c(alpha.quant, -alpha.quant)
    CI.cmle <- cf[2,1] - cf[2,2] * c(alpha.quant, -alpha.quant)
    CI.rp <- predict(smooth.spline(r.p$y, r.p$x),
                     c(alpha.quant, -alpha.quant), 0)$y
    CI.rmp <- predict(smooth.spline(r.mp$y, r.mp$x),
                      c(alpha.quant, -alpha.quant), 0)$y
    CI.lr.mp <- predict(smooth.spline(lr.mp$y, lr.mp$x),
                        c(1-alpha/2, 1-alpha, alpha/2, alpha), 0)$y
    if(is.scale)
    {
      CI.mle <- exp(CI.mle)
      CI.cmle <- exp(CI.cmle)
      CI.rp <- exp(CI.rp)
      CI.rmp <- exp(CI.rmp)
      CI.lr.mp <- exp(CI.lr.mp)
    }
    CI <- t(matrix(c(CI.mle, CI.cmle, CI.rp, CI.rmp, CI.lr.mp),
                   ncol = 4*dim.alpha, byrow = TRUE))
    dimnames(CI) <- list(c(paste("Lower conf. bound (2sd,", 100*(1-alpha), "%)"),
                           paste("Lower conf. bound (1sd,", 100*(1-alpha), "%)"),
                           paste("Upper conf. bound (2sd,", 100*(1-alpha), "%)"),
                           paste("Upper conf. bound (1sd,", 100*(1-alpha), "%)")),
                         c("Wald pivot.              ",
                           "Wald pivot (cond.)       ",
                           "Likelihood root          ",
                           "Modified likelihood root ",
                           "Lugannani-Rice approx.   ") )
  }
  old.test <- test
  if(!is.null(test))
  {
    if(is.scale) test <- log(test)
    test.mle <- (cf[1,1] - test) / cf[1,2]
    test.cmle <- (cf[2,1] - test) / cf[2,2]
    test.rp <- predict(smooth.spline(r.p), test, 0)$y
    test.rmp <- predict(smooth.spline(r.mp), test, 0)$y
    test.lr.mp <- predict(smooth.spline(lr.mp), test, 0)$y
    test.lr.mp <- ifelse(test.lr.mp > 0.5, 1 - test.lr.mp, test.lr.mp)
    TEST <- t(matrix(c(test.mle,
                       ifelse(test.mle > 0, (1 - pnorm(test.mle)),
                                            pnorm(test.mle)),
                       test.cmle,
                       ifelse(test.cmle > 0, (1-pnorm(test.cmle)),
                                             pnorm(test.cmle)),
                       test.rp,
                       ifelse(test.rp > 0, (1-pnorm(test.rp)),
                                           pnorm(test.rp)),
                       test.rmp,
                       ifelse(test.rmp > 0, (1-pnorm(test.rmp)),
                                            pnorm(test.rmp)),
                       rep(NA, dim.test), test.lr.mp
                      ),
                      nrow = 2 * dim.test))
     dimnames(TEST) <- list(
                         c("Wald pivot               ",
                           "Wald pivot (cond.)       ",
                           "Likelihood root          ",
                           "Modified likelihood root ",
                           "Lugannani-Rice approx.   "),
              c(paste("statistic (", object$offset, "=", test, ")  "),
                paste("tail prob. (", object$offset, "=", test, 
                                                         ")        ")))
     qT <- object$workspace$q
     qT <- predict(smooth.spline(qT, all.knots=FALSE), test, 0)$y
     qT <- matrix(qT, ncol=1)
     dimnames(qT) <- list(paste("(", object$offset, "=", test, ")"),
                          "q")
    TEST <- list(stats=TEST, qTerm=qT)
  }
  summary.object <- list( coefficients = object$coefficients,
                          conf.int = if(int) CI,
                          signif.tests = if(!is.null(test)) TEST ,
                          call = object$call,
                          formula = object$formula,
                          family = object$family,
                          offset = object$offset,
                          alpha = alpha,
                          hypotheses = old.test,
                          diagnostics = object$diagnostics,
                          n.approx = object$n.approx,
                          all = all, cf = coef, int = int,
                          is.scalar = is.scalar, digits = digits )
  class(summary.object) <- c("summary.marg", "summary.cond")
  summary.object
}

print.summary.marg <- function(x, all = x$all, Coef = x$cf,
                               int = x$int, test = x$hyp,
                               digits = if(!is.null(x$digits)) x$digits
                                        else max(3, getOption("digits")-3),
                               ...)
{
  cat("\n Formula:  ")
  dput(x$formula)
  cat(" Family:  ")
  print(x$family)
  cat(" Offset:  ")
  if( x$offset != "1" )
    print(as.name(x$offset))
  else  cat("Intercept\n")
  cat("\n")
  is.scale <- ( paste(x$offset) == "scale" )
  is.scalar <- x$is.scalar
  if(is.scalar)
    cat("Original model contains no nuisance parameter\n\n")
  if(Coef)
  {
    if( !is.scalar )
      print(signif(x$coefficients, digits=digits))
    else  print(signif(x$coefficients[1,], digits=digits))
  }
  if(int)
  {
    dim.alpha <- length(x$alpha)
    names.stat <- c("Wald pivot               ",
                    "Wald pivot (cond.)       ",
                    "Likelihood root          ",
                    "Modified likelihood root ",
                if(all)
                  c("Lugannani-Rice approx.   "))
    cat("\nConfidence intervals",
        if(all) c("and one-sided confidence bounds") )
    cat("\n--------------------",
        if(all) c("-------------------------------"),
        sep="-")

    for(i in seq(along=x$alpha))
    {
      cat("\n level =", 100*(1-x$alpha[i]), "%\n")
      if(all)
      {
         n <- length(names.stat)
         idx <- 1:n
         ntimes <- n
         if( is.scalar )
         {
           idx <- idx[-2]
           ntimes <- ntimes - 1
         }
         print.mat <- data.frame(
                 as.vector(signif(x$conf.int[i,idx],
                                  digits=digits)), rep("", ntimes),
                 as.vector(signif(x$conf.int[i+2*dim.alpha ,idx],
                                  digits=digits)))
         dimnames(print.mat)[[2]] <-
            c("lower", "two-sided", "upper")
         dimnames(print.mat)[[1]] <- names.stat[idx]
         print(print.mat)
         cat("\n")
         print.mat <- data.frame(
                 as.vector(signif(x$conf.int[dim.alpha+i,idx],
                                  digits=digits)), rep("", ntimes),
                 as.vector(signif(x$conf.int[i+3*dim.alpha ,idx],
                                  digits=digits)))
         dimnames(print.mat)[[2]] <-
            c("lower", " one-sided", "upper")
         dimnames(print.mat)[[1]] <- names.stat[idx]
         print(print.mat)
      }
      else
      {
         n <- length(names.stat)
         idx <- 1:n
         ntimes <- n
         if( is.scalar )
         {
           idx <- idx[-2]
           ntimes <- ntimes - 1
         }
         print.mat <- data.frame(
                 as.vector(signif(x$conf.int[i,idx], 
                                  digits=digits)), rep("", ntimes),
                 as.vector(signif(x$conf.int[i+2*dim.alpha ,idx],
                                  digits=digits)))
         dimnames(print.mat)[[2]] <-
            c("lower", "two-sided", "upper")
         dimnames(print.mat)[[1]] <- names.stat[idx]
         print(print.mat)
      }
    }
  }
  qT <- x$signif.tests$qTerm
  if( !is.null(test) )
  {
    dim.test <- length(x$hypotheses)
    names.stat <- c("Wald pivot               ",
                    "Wald pivot (cond.)       ",
                    "Likelihood root          ",
                    "Modified likelihood root ",
                if(all)
                  c("Lugannani-Rice approx.   ") )
    cat("\nTest statistics")
    cat("\n---------------")
    for(i in seq(along = x$hypotheses))
    {
      mat <- matrix(signif(x$signif.tests$stats[1:(ifelse(all, 5, 4)),
                                               c(i, i+dim.test)],
                           digits = digits), ncol = 2)
      dimnames(mat) <- list(names.stat, c(" statistic",
                                          " tail prob."))
      txt <- if(is.scale) "\n hypothesis : scale =" 
             else if(x$offset == "1")
                    "\n hypothesis : Intercept ="
                  else paste("\n hypothesis : coef(", x$offset, ") =")
      txt <- paste(txt, x$hypotheses[i], "\n")
      cat(txt)
      if( !is.scalar )
        print(mat)
      else  print(mat[-2,])
      xq <- matrix(qT[i,], nrow=1)
      dimnames(xq) <- list("\"q\" correction term:", "")
      print(xq, digits=digits)
    } 
  }
  xd <- x$diagnostics 
  cat("\nDiagnostics:")
  cat("\n----------- \n")
  if( !is.scalar )
    print(xd, digits=digits)
  else  print(xd[1], digits=digits)
  cat("\n Approximation based on", x$n.approx, "points\n")
}
