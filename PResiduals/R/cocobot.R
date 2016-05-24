#'  Conditional continuous by ordinal tests for association.
#' 
#' \code{cocobot} tests for independence between an ordered categorical
#' variable, \var{X}, and a continuous variable, \var{Y}, conditional on other variables,
#' \var{Z}.  The basic approach involves fitting an ordinal model of \var{X} on
#' \var{Z}, a linear model of \var{Y} on \var{Z}, and then determining whether there is any
#' residual information between \var{X} and \var{Y}.  This is done by
#' computing residuals for both models, calculating their correlation, and 
#' testing the null of no residual correlation.  This procedure is analogous to test statistic 
#' \code{T2} in \code{cobot}.  Two test statistics (correlations) are currently output.  The first
#' is the correlation between probability-scale residuals. The second is the correlation between 
#' the observed-minus-expected residual for the continuous outcome model and a latent variable residual
#' for the ordinal model.
#' 
#' Formula is specified as \code{\var{X} | \var{Y} ~ \var{Z}}.
#' This indicates that models of \code{\var{X} ~ \var{Z}} and
#' \code{\var{Y} ~ \var{Z}} will be fit.  The null hypothsis to be
#' tested is \eqn{H_0 : X}{H0 : X} independant of \var{Y} conditional
#' on \var{Z}.  The ordinal variable, \code{\var{X}}, must precede the \code{|} and be a factor variable, and \code{\var{Y}} must be continuous.
#'  
#' @references Li C and Shepherd BE (2012) 
#' A new residual for ordinal outcomes.
#' \emph{Biometrika}. \bold{99}: 473--480.
#' @references Shepherd BE, Li C, Liu Q (submitted)
#' Probability-scale residuals for continuous, discrete, and censored data.
#'

#' @param formula an object of class \code{\link{Formula}} (or one
#' that can be coerced to that class): a symbolic description of the
#' model to be fitted.  The details of model specification are given
#' under \sQuote{Details}.
#' 
#' @param link The link family to be used for the ordinal model of 
#' \var{X} on \var{Z}.  Defaults to \samp{logit}. Other options are
#' \samp{probit}, \samp{cloglog}, and \samp{cauchit}.
#' 
#' @param data an optional data frame, list or environment (or object
#' coercible by \code{\link{as.data.frame}} to a data frame)
#' containing the variables in the model.  If not found in
#' \code{data}, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{cocobot} is called.
#' 
#' @param subset an optional vector specifying a subset of
#' observations to be used in the fitting process.
#' 
#' @param na.action action to take when \code{NA} present in data.
#' 
#' @param emp logical indicating whether the residuals from the model of
#' \var{Y} on \var{Z} are computed based on the assumption of normality (\code{FALSE}) 
#' or empirically (\code{TRUE}).
#' 
#' @param fisher logical indicating whether to apply fisher transformation to compute confidence intervals and p-values for the correlation.
#' 
#' @param conf.int numeric specifying confidence interval coverage.
#' 
#' @return object of \samp{cocobot} class.
#' @export
#' @examples
#' data(PResidData)
#' cocobot(y|w ~ z, data=PResidData)
#' @importFrom stats qlogis qnorm qcauchy integrate

cocobot <- function(formula, data, link=c("logit", "probit", "cloglog", "cauchit"),
                      subset, na.action=getOption('na.action'), 
                      emp=TRUE,fisher=FALSE,conf.int=0.95) {

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
  
  score.xz <- ordinal.scores(mx, mxz, method=link)
  if(emp == TRUE) {
      score.yz <- lm.scores(y=model.response(my), X=myz, emp=TRUE)
  }
  else {
      score.yz <- lm.scores(y=model.response(my), X=myz, emp=FALSE)
  }
  
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

  ### presid vs obs-exp
  #ta <-  corTS(xz.presid, score.yz$resid,
  #              score.xz$dl.dtheta, score.yz$dl.dtheta,
  #              score.xz$d2l.dtheta.dtheta, score.yz$d2l.dtheta.dtheta,
  #              xz.dpresid.dtheta, score.yz$dresid.dtheta, fisher)
  #ans$TS$TA <- 
  #    list( 
  #      ts=ta$TS, var=ta$var.TS, pval=ta$pval.TS,
  #      label='PResid vs. Obs-Exp'
  #    )
  
  ## presid vs presid (use pdf of normal)
  tb <- corTS(xz.presid, score.yz$presid,
              score.xz$dl.dtheta, score.yz$dl.dtheta,
              score.xz$d2l.dtheta.dtheta, score.yz$d2l.dtheta.dtheta,
              xz.dpresid.dtheta, score.yz$dpresid.dtheta,fisher)
  if (emp==TRUE){
    tb.label = "PResid vs. PResid (empirical)"
  } else {
      tb.label = "PResid vs. PResid (assume normality)"    
  }
  ans$TS$TB <- 
      list( 
        ts=tb$TS, var=tb$var.TS, pval=tb$pval.TS,
        label = tb.label
      )

  T3 <- 3*sum(xz.presid*score.yz$presid)/N
  
  rij <- cbind(score.xz$Gamma, 1)[cbind(1:N, xx)]
  rij_1 <- cbind(0,score.xz$Gamma)[cbind(1:N, xx)]
  pij <- rij-rij_1

  G.inverse <- switch(link[1], logit = qlogis, probit = qnorm,
                 cloglog = qgumbel, cauchit = qcauchy)
  xz.latent.resid <- rep(NA, N)
  
  inverse_fail <- FALSE 
  for (i in 1:N){
    tmp <- try(integrate(G.inverse, rij_1[i], rij[i])$value/pij[i],silent=TRUE)
    if (inherits(tmp,'try-error')){
      if (link[1] != 'cauchit')
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
    
    
    ### latent.resid vs obs-exp
    tc <- corTS(xz.latent.resid, score.yz$resid,
                 score.xz$dl.dtheta, score.yz$dl.dtheta,
                 score.xz$d2l.dtheta.dtheta, score.yz$d2l.dtheta.dtheta,
                 xz.dlatent.dtheta, score.yz$dresid.dtheta,fisher)
    ans$TS$TC <- 
        list( 
          ts=tc$TS, var=tc$var.TS, pval=tc$pval.TS,
          label = 'Latent.resid vs. Obs-Exp'
        )
    
    #if (emp==TRUE){
    #  ### latent vs presid (emprical)
    #  td = corTS(xz.latent.resid, score.yz$presid.k,
    #              score.xz$dl.dtheta, score.yz$dl.dtheta,
    #              score.xz$d2l.dtheta.dtheta, score.yz$d2l.dtheta.dtheta,
    #              xz.dlatent.dtheta, score.yz$dpresid.dtheta.k,fisher)
    #  td.label = "Latent.resid  vs. PResid (empirical)"
    #} else {
    #  ### latent vs presid (use pdf of normal)
    #  td = corTS(xz.latent.resid, score.yz$presid,
    #              score.xz$dl.dtheta, score.yz$dl.dtheta,
    #              score.xz$d2l.dtheta.dtheta, score.yz$d2l.dtheta.dtheta,
    #              xz.dlatent.dtheta, score.yz$dpresid.dtheta,fisher)
    #  td.label = "Latent.resid  vs. PResid (assume normality)"
    #  
    #}
    #ans$TS$TD <- 
    #    list( 
    #      ts=td$TS, var=td$var.TS, pval=td$pval.TS,
    #      label = td.label
    #    )
  }
  
  ans <- structure(ans, class="cocobot")

  # Apply confidence intervals
  for (i in seq_len(length(ans$TS))) {
    ts_ci <- getCI(ans$TS[[i]]$ts,ans$TS[[i]]$var,ans$fisher,conf.int)
    ans$TS[[i]]$lower <- ts_ci[1]
    ans$TS[[i]]$upper <- ts_ci[2]
  }

  print(T3)
  return(ans)
}
