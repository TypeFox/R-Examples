#'  Conditional tests for association.
#' 
#' \code{megabot} tests for correlation between a variable, \var{X}, and another variable, \var{Y},
#' conditional on other variables, \var{Z}.  The basic approach involves fitting an specified model of \var{X} on
#' \var{Z}, a specified model of \var{Y} on \var{Z}, and then determining whether there is any
#' residual information between \var{X} and \var{Y}.  This is done by
#' computing residuals for both models, calculating their correlation, and 
#' testing the null of no residual correlation.  The test statistic output 
#' is the correlation between probability-scale residuals.  \var{X} and \var{Y} can
#' be continous or ordered discrete variables.  \code{megabot} replicates the functionality
#' of \code{\link{cobot}}, \code{\link{cocobot}}, and \code{\link{countbot}}
#' 
#' Formula is specified as \code{\var{X} | \var{Y} ~ \var{Z}}.
#' This indicates that models of \code{\var{X} ~ \var{Z}} and
#' \code{\var{Y} ~ \var{Z}} will be fit.  The null hypothesis to be
#' tested is \eqn{H_0 : X}{H0 : X} independent of \var{Y} conditional
#' on \var{Z}.
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
#' @param fit.x,fit.y The fitting function used for the model of \var{X} or \var{Y} on
#' \var{Z}.  Options are \samp{ordinal}, \samp{lm}, \samp{lm.emp},
#' \samp{poisson}, \samp{nb}, and \samp{orm}.
#' 
#' @param link.x,link.y The link family to be used for the ordinal model of 
#' \var{X} on \var{Z}.  Defaults to \samp{logit}. Other options are
#' \samp{probit}, \samp{cloglog}, and \samp{cauchit}. Used only when
#' \samp{fit.x} is either \samp{ordinal} or \samp{orm}.
#' 
#' @param data an optional data frame, list or environment (or object
#' coercible by \code{\link{as.data.frame}} to a data frame)
#' containing the variables in the model.  If not found in
#' \code{data}, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{megabot} is called.
#' 
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
#' @examples
#' data(PResidData)
#' megabot(y|w ~ z, fit.x="ordinal", fit.y="lm.emp", data=PResidData)
#' @importFrom rms lrm


megabot <- function(formula, data, fit.x, fit.y, link.x=c("logit", "probit", "cloglog", "cauchit"),
                    link.y = c("logit", "probit", "cloglog", "cauchit"),
                    subset, na.action=getOption('na.action'), fisher=FALSE,conf.int=0.95){

    ## Construct the model frames for x ~ z and y ~ z
    F1 <- Formula(formula)
    Fx <- formula(F1, lhs=1)
    Fy <- formula(F1, lhs=2)
    
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
                 "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf$na.action <- na.action
    ## We set xlev to a benign non-value in the call so that it won't get partially matched
    ## to any variable in the formula. For instance a variable named 'x' could possibly get
    ## bound to xlev, which is not what we want.
    mf$xlev <- integer(0) 
    mf[[1L]] <- as.name("model.frame")
  
  
    mx <- my <- mf

    ## NOTE: we add the opposite variable to each model frame call so that
    ## subsetting occurs correctly. Later we strip them off.
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

    ## Construct the model matrix z
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

    score.xz <- switch(fit.x,
                       ordinal = ordinal.scores(mx, mxz, method=link.x),
                       lm = lm.scores(y=model.response(mx), X=mxz, emp=FALSE),
                       lm.emp = lm.scores(y=model.response(mx), X=mxz, emp=TRUE),
                       poisson = poisson.scores(y=model.response(mx), X=mxz),
                       nb = nb.scores(y=model.response(mx), X=mxz),
                       orm = orm.scores(y=model.response(mx), X=mxz, link=link.x))
    
    score.yz <- switch(fit.y,
                       ordinal = ordinal.scores(my, myz, method=link.y),
                       lm = lm.scores(y=model.response(my), X=myz, emp=FALSE),
                       lm.emp = lm.scores(y=model.response(my), X=myz, emp=TRUE),
                       poisson = poisson.scores(y=model.response(my), X=myz),
                       nb = nb.scores(y=model.response(my), X=myz),
                       orm = orm.scores(y=model.response(my), X=myz, link=link.y))


    ## return value
    ans <- list(
        TS=list(),
        fisher=fisher,
        conf.int=conf.int,
        data.points=data.points
        )

    ## presid vs presid (use pdf of normal)
    tb = corTS(score.xz$presid, score.yz$presid,
        score.xz$dl.dtheta, score.yz$dl.dtheta,
        as.matrix(score.xz$d2l.dtheta.dtheta), as.matrix(score.yz$d2l.dtheta.dtheta),
        score.xz$dpresid.dtheta, score.yz$dpresid.dtheta,fisher)
    tb.label = "PResid vs. PResid (assume normality)"
    

    ans$TS$TB <- list( 
        ts=tb$TS, var=tb$var.TS, pval=tb$pval.TS,
        label = tb.label
        )

    ans <- structure(ans, class="cocobot")

    ## Apply confidence intervals
    for (i in seq_len(length(ans$TS))) {
        ts_ci <- getCI(ans$TS[[i]]$ts,ans$TS[[i]]$var,ans$fisher,conf.int)
        ans$TS[[i]]$lower <- ts_ci[1]
        ans$TS[[i]]$upper <- ts_ci[2]
    }

    return(ans)
}
