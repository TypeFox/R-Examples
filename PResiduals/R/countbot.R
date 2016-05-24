#' Conditional count by ordinal tests for association.
#'
#' \code{countbot} tests for independence between an ordered categorical
#' variable, \var{X}, and a count variable, \var{Y}, conditional on other variables,
#' \var{Z}.  The basic approach involves fitting an ordinal model of \var{X} on
#' \var{Z}, a Poisson or Negative Binomial model of \var{Y} on \var{Z}, and then determining whether there is any
#' residual information between \var{X} and \var{Y}.  This is done by
#' computing residuals for both models, calculating their correlation, and 
#' testing the null of no residual correlation.  This procedure is analogous to test statistic 
#' \code{T2} in \code{cobot}.  Two test statistics (correlations) are currently output.  The first
#' is the correlation between probability-scale residuals. The second is the correlation between 
#' the Pearson residual for the count outcome model and a latent variable residual
#' for the ordinal model.

#'
#' Formula is specified as \code{\var{X} | \var{Y} ~ \var{Z}}.
#' This indicates that models of \code{\var{X} ~ \var{Z}} and
#' \code{\var{Y} ~ \var{Z}} will be fit.  The null hypothesis to be
#' tested is \eqn{H_0 : X}{H0 : X} independent of \var{Y} conditional
#' on \var{Z}.  The ordinal variable, \code{\var{X}}, must precede the \code{|} and be a factor variable, and \code{\var{Y}} must be an integer.
#' @references Li C and Shepherd BE (2012) 
#' A new residual for ordinal outcomes.
#' \emph{Biometrika}. \bold{99}: 473--480.
#' @references Shepherd BE, Li C, Liu Q (submitted)
#' Probability-scale residuals for continuous, discrete, and censored data.
#'

#'
#' @param formula an object of class \code{\link{Formula}} (or one
#' that can be coerced to that class): a symbolic description of the
#' model to be fitted.  The details of model specification are given
#' under \sQuote{Details}.
#'
#' @param data an optional data frame, list or environment (or object
#' coercible by \code{\link{as.data.frame}} to a data frame)
#' containing the variables in the model.  If not found in
#' \code{data}, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{countbot} is called.

#' @param link.x The link family to be used for the ordinal model of 
#' \var{X} on \var{Z}.  Defaults to \samp{logit}. Other options are
#' \samp{probit}, \samp{cloglog}, and \samp{cauchit}.
#' 
#' @param family.y The error distribution for the count model of \var{Y} on \var{Z}.
#' Defaults to \samp{poisson}. The other option is \samp{negative binomial}. 
#' If \samp{negative binomial} is specified, \code{\link[MASS]{glm.nb}} is called to fit the count model.

#' @param subset an optional vector specifying a subset of
#' observations to be used in the fitting process.
#' 
#' @param na.action action to take when \code{NA} present in data.
#' 
#' @param fisher logical indicating whether to apply fisher transformation to compute confidence intervals and p-values for the correlation.
#' 
#' @param conf.int numeric specifying confidence interval coverage.
#' 
#' @return object of \samp{cocobot} class.
#' @export
#' @importFrom stats qlogis qnorm qcauchy integrate
#' @examples
#'
#' generate.data3 = function(alphax, betax, alphay, betay, eta, N) {
#'   z = rnorm(N,0,1)
#'   x = v = numeric(N)
#'  
#'   ## px is an N x length(alphax) matrix.
#'   ## Each row has the TRUE cummulative probabilities for each subject.
#'   px = (1 + exp(- outer(alphax, betax*z, "+"))) ^ (-1)
#'   aa = runif(N)
#'   for(i in 1:N)
#'     x[i] = sum(aa[i] > px[,i])
#'   x = as.numeric(as.factor(x))
#'   v = rpois(N, exp(outer(alphay, betay*z+eta[x], "+")))
#'  
#'   return(list(x=as.factor(x), v=v, z=z))
#' }
#'
#' set.seed(13)
#' alphax = c(-1, 0, 1, 2)
#' betax = 1
#' alphay = 1
#' betay = -.5
#'
#' #eta = rep(0, 5)
#' eta = c(1:5)/5
#' N = 100
#'
#' data <- generate.data3(alphax, betax, alphay, betay, eta, N)
#' ### 
#' countbot(x|v~z, data=data, link.x="logit", family.y="poisson")
#' countbot(x|v~z, data=data, link.x="probit", family.y="negative binomial")





countbot <- function(formula, data, link.x=c("logit", "probit", "cloglog", "cauchit"),
                     family.y=c("poisson", "negative binomial"),
                     subset, na.action=getOption('na.action'), 
                     fisher=FALSE, conf.int=0.95) {
  
  
  # Construct the model frames for x ~ z and y ~ z
  F1 <- Formula(formula)
  Fx <- formula(F1, lhs=1)
  Fy <- formula(F1, lhs=2)
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.action
  # We set xlev to a benign non-value in the call so that it won't get partially matched
  # to any variable in the formula. For instance a variable named 'x' could possibly get
  # bound to xlev, which is not what we want.
  mf$xlev <- integer(0) 
  mf[[1L]] <- as.name("model.frame")
  
  
  mx <- my <- mf
  
  # NOTE: we add the opposite variable to each model frame call so that
  # subsetting occurs correctly. Later we strip them off.
  mx[["formula"]] <- Fx
  yName <- all.vars(Fy[[2]])[1]
  mx[[yName]] <- Fy[[2]]
  
  my[["formula"]] <- Fy
  xName <- all.vars(Fx[[2]])[1]
  my[[xName]] <- Fx[[2]]
  
  mx <- eval(mx, parent.frame())
  mx[[paste('(',yName,')',sep='')]] <- NULL
  
  my <- eval(my, parent.frame())
  my[[paste('(',xName,')',sep='')]] <- NULL
  
  data.points <- nrow(mx)
  
  if (!is.factor(mx[[1]])){
    warning("Coercing ",names(mx)[1]," to factor. Check the ordering of categories.")
    mx[[1]] <- as.factor(mx[[1]])
  }
  
  if (is.factor(my[[1]])){
    stop(names(my)[1]," cannot be a factor.")
  }
  
  # Construct the model matrix z
  mxz <- model.matrix(attr(mx,'terms'),mx) 
  zzint <- match("(Intercept)", colnames(mxz), nomatch = 0L)
  if(zzint > 0L) {
    mxz <- mxz[, -zzint, drop = FALSE]
  }
  
  myz <- model.matrix(attr(my,'terms'),my) 
  zzint <- match("(Intercept)", colnames(myz), nomatch = 0L)
  if(zzint > 0L) {
    myz <- myz[, -zzint, drop = FALSE]
  }
  
  score.xz <- ordinal.scores(mx, mxz,method=link.x[1])
  if (family.y[1]=="poisson")
    score.yz <- poisson.scores(y=model.response(my), X=myz)
  else if (family.y[1]=="negative binomial")
    score.yz <- nb.scores(y=model.response(my), X=myz)
  else stop("family.y has to be 'poisson' or 'negative binomial'")
  
  npar.xz = dim(score.xz$dl.dtheta)[2]
  npar.yz = dim(score.yz$dl.dtheta)[2]
  xx = as.integer(model.response(mx))
  
  nx = length(table(xx))
  
  N = length(xx)
  low.x = cbind(0, score.xz$Gamma)[cbind(1:N, xx)]
  hi.x = cbind(1-score.xz$Gamma, 0)[cbind(1:N, xx)]
  
  xz.presid <- low.x - hi.x
  xz.dpresid.dtheta <- score.xz$dlow.dtheta - score.xz$dhi.dtheta
  
  ## return value
  ans <- list(
    TS=list(),
    fisher=fisher,
    conf.int=conf.int,
    data.points=data.points
  )
  
  tb = corTS(xz.presid, score.yz$presid,
             score.xz$dl.dtheta, score.yz$dl.dtheta,
             score.xz$d2l.dtheta.dtheta, score.yz$d2l.dtheta.dtheta,
             xz.dpresid.dtheta, score.yz$dpresid.dtheta,fisher)
  tb.label = "PResid vs. PResid" 
  
  ans$TS$TB <- 
    list( 
      ts=tb$TS, var=tb$var.TS, pval=tb$pval.TS,
      label = tb.label
    )

  T3 <- 3 * sum(xz.presid * score.yz$presid) / N
  rij <- cbind(score.xz$Gamma, 1)[cbind(1:N, xx)]
  rij_1 <- cbind(0,score.xz$Gamma)[cbind(1:N, xx)]
  pij <- rij-rij_1
  
  G.inverse <- switch(link.x[1], logit = qlogis, probit = qnorm,
                      cloglog = qgumbel, cauchit = qcauchy)
  xz.latent.resid <- rep(NA, N)
  
  inverse_fail <- FALSE 
  for (i in 1:N){
    tmp <- try(integrate(G.inverse, rij_1[i], rij[i])$value/pij[i],silent=TRUE)
    if (inherits(tmp,'try-error')){
      if (link.x[1] != 'cauchit')
        warning("Cannot compute latent variable residual.")
      else
        warning("Cannot compute latent variable residual with link function cauchit.")
      inverse_fail <- TRUE
      break
    } else {
      xz.latent.resid[i] <- tmp
    }
  }
  
  if (!inverse_fail){
    ### To compute dlatent.dtheta (need dgamma.dtheta and dp0.dtheta from ordinal scores)
    xz.dlatent.dtheta = dpij.dtheta = matrix(, npar.xz, N)
    
    drij_1.dtheta <- score.xz$dlow.dtheta
    drij.dtheta <- -score.xz$dhi.dtheta
    for(i in 1:N) {  
      dpij.dtheta[,i] <- score.xz$dp0.dtheta[i, xx[i],]
      
      if (xx[i] == 1) {
        xz.dlatent.dtheta[,i] <- -xz.latent.resid[i]/pij[i]*dpij.dtheta[,i] + 1/pij[i]*(
          G.inverse(rij[i])*drij.dtheta[,i] - 0 )
      } else if(xx[i] == nx){
        xz.dlatent.dtheta[,i] <- -xz.latent.resid[i]/pij[i]*dpij.dtheta[,i] + 1/pij[i]*(
          0 - G.inverse(rij_1[i])*drij_1.dtheta[,i] )
      } else
        xz.dlatent.dtheta[,i] <- -xz.latent.resid[i]/pij[i]*dpij.dtheta[,i] + 1/pij[i]*(
          G.inverse(rij[i])*drij.dtheta[,i] - G.inverse(rij_1[i])*drij_1.dtheta[,i])
    }
    
    
    ### latent.resid vs pearson resid
    tc <- corTS(xz.latent.resid, score.yz$pearson.resid,
                score.xz$dl.dtheta, score.yz$dl.dtheta,
                score.xz$d2l.dtheta.dtheta, score.yz$d2l.dtheta.dtheta,
                xz.dlatent.dtheta, score.yz$dpearson.resid.dtheta, fisher)
    
    ans$TS$TC <- 
      list( 
        ts=tc$TS, var=tc$var.TS, pval=tc$pval.TS,
        label = 'Latent.resid vs. Pearson.resid'
      )
  }
  ans <- structure(ans, class="cocobot")
  
  # Apply confidence intervals
  for (i in seq_len(length(ans$TS))){
    ts_ci <- getCI(ans$TS[[i]]$ts,ans$TS[[i]]$var,ans$fisher,conf.int)
    ans$TS[[i]]$lower <- ts_ci[1]
    ans$TS[[i]]$upper <- ts_ci[2]
  }

  print(T3)
  return(ans)
  
}


#### example
## generate count by ordinal data
## generate.data3 = function(alphax, betax, alphay, betay, eta, N) {
##   z = rnorm(N,0,1)
##   x = y = numeric(N)
  
##   ## px is an N x length(alphax) matrix.
##   ## Each row has the TRUE cummulative probabilities for  each subject.
##   px = (1 + exp(- outer(alphax, betax*z, "+"))) ^ (-1)
##   aa = runif(N)
##   for(i in 1:N)
##     x[i] = sum(aa[i] > px[,i])
##   x = as.numeric(as.factor(x))
##   y = rpois(N, exp(outer(alphay, betay*z+eta[x], "+")))
  
##   return(list(x=as.factor(x), y=y, z=z))
## }


## set.seed(13)
## alphax = c(-1, 0, 1, 2)
## betax = 1
## alphay = 1
## betay = -.5

## #eta = rep(0, 5)
## eta = c(1:5)/20 
## N = 200
## data <- generate.data3(alphax, betax, alphay, betay, eta, N)

## #### check for cocobot
## cocobot(x|y~z, data=data)
## countbot(x|y~z, data=data, fisher=TRUE)
## countbot(x|y~z, data=data, family="negative binomial")
