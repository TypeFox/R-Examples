#' Print continuous ordinal regression objects
#'
#' \code{print.ocm} is the ordinalCont specific method for the generic function \code{print}, 
#' which prints objects of class \code{ocm}.
#' @param x an object of class \code{ocm}, usually, a result of a call to \code{ocm}
#' @param ... further arguments passed to or from other methods
#' @return Prints an \code{ocm} object
#' @keywords likelihood, log-likelihood.
#' @method print ocm
#' @seealso \code{\link{ocm}}, \code{\link{summary.ocm}}
#' @examples 
#' ANZ0001.ocm <- ANZ0001[ANZ0001$cycleno==0 | ANZ0001$cycleno==5,]
#' ANZ0001.ocm$cycleno[ANZ0001.ocm$cycleno==5] <- 1
#' fit.overall  <- ocm(overall  ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' print(fit.overall)
#' @export
#' @author Maurizio Manuguerra, Gillian Heller

print.ocm <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients, ...)
}

#' @title Summarizing Continuous Ordinal Fits
#' @description Summary method for class \code{ocm}
#' @param object an object of class \code{ocm}, usually a result of a call to \code{ocm}
#' @param ... further arguments passed to or from other methods
#' @method summary ocm
#' @keywords summary
#' @seealso \code{\link{ocm}}, \code{\link{print.ocm}}
#'  @author Maurizio Manuguerra, Gillian Heller
#' @examples
#' ANZ0001.ocm <- ANZ0001[ANZ0001$cycleno==0 | ANZ0001$cycleno==5,]
#' ANZ0001.ocm$cycleno[ANZ0001.ocm$cycleno==5] <- 1
#' fit.overall  <- ocm(overall  ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' summary(fit.overall)
#' @export

summary.ocm <- function(object, ...)
{
  se <- sqrt(diag(object$vcov))
  tval <- coef(object)[1:length(se)] / se
  TAB <- data.frame(Estimate = coef(object),
               StdErr = se,
               t.value = tval,
               p.value = 2*pt(-abs(tval), df=object$df))
  res <- list(call=object$call,
              coefficients=TAB,
              len_beta=object$len_beta,
              len_gfun=object$len_gfun)
  class(res) <- "summary.ocm"
  print(res, ...)
}


print.summary.ocm <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Coefficients:\n")
  printCoefmat(x$coefficients[1:x$len_beta,], P.values = TRUE, has.Pvalue = TRUE, signif.legend = FALSE, ...)
  cat("\n")
  cat("g function:\n")
  printCoefmat(x$coefficients[(x$len_beta+1):(x$len_beta+x$len_gfun),], P.values = TRUE, has.Pvalue = TRUE, ...)
}




#' @title Predict method for Continuous Ordinal Fits
#' 
#' @description Predicted values based on \code{ocm} object
#' @param object an object of class \code{ocm}, usually a result of a call to \code{ocm}
#' @param newdata optionally, a data frame in which to look for variables with 
#' which to predict. 
#' Note that all predictor variables should be present, having the same names as the variables 
#' used to fit the model. If \code{NULL}, predictions are computed for the original dataset.
#' @param ndens the number of points on the continuous ordinal scale (0, 1) over which the densities are computed. 
#' The default is 100.
#' @param ... further arguments passed to or from other methods
#' @keywords predict
#' @method predict ocm
#' @author Maurizio Manuguerra, Gillian Heller
#' @return  A list containing the following components: 
#' \item{mode}{a vector of length equal to the number of observations.
#' Each element is the mode of \code{v}, 
#' the  continuous ordinal random variable, conditional on the covariates in the model.}
#' \item{density}{a matrix with number of rows equal to the number of observations. Each row 
#' contains the values of the density function of \code{v} conditional on the covariates in the 
#' model. 
#' The density function is calculated over \code{ndens} equally-spaced values of v in (0,1).}
#' \item{x}{a vector with the \code{ndens} equally-spaced values of \code{v} in (0,1) used to compute the 
#' density of v}
#' \item{formula}{the formula used to fit the model}
#' \item{newdata}{a new data frame used to make predictions. It takes value NULL if no new data frame has been used.}
#' @details An object of class \code{ocm} and optionally a new data 
#' frame are used to compute the probability 
#' densities of \code{v}, the continuous ordinal score. The estimated parameters 
#' of the fitted model and \code{ndens} (default: 100) 
#' values of \code{v} are used to compute the probability densities on the latent scale. 
#' These values are then transformed to scores on the continuous ordinal 
#' scale using the g function and the estimated values 
#' of \code{M}, \code{B}, and \code{T}.
#' @examples 
#' ANZ0001.ocm <- ANZ0001[ANZ0001$cycleno==0 | ANZ0001$cycleno==5,]
#' ANZ0001.ocm$cycleno[ANZ0001.ocm$cycleno==5] <- 1
#' fit.overall <- ocm(overall ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' pred <- predict(fit.overall)
#' plot(pred)
#' @seealso \code{\link{ocm}}, \code{\link{plot.predict.ocm}}
#' @export


predict.ocm <- function(object, newdata=NULL, ndens=100, ...)
{
  formula <- object$formula
  params <- coef(object)
  if(is.null(newdata)){
    x <- object$x 
  }else{
    x <- model.matrix(object$formula, newdata)
  }
  len_beta <- ncol(x)
  v <- seq(0, 1, length.out = ndens+2)[2:(ndens+1)]
  modes <- NULL
  densities <- NULL
  #FIXME: rewrite efficiently
  for (subject in 1:nrow(x)){
    d.matrix <- matrix(rep(x[subject,], ndens), nrow = ndens, dimnames = list(as.character(1:ndens), 
                                                                              colnames(x)), byrow = TRUE)
    densities <- rbind(densities, t(logdensity_glf(par = params, v = v, d.matrix = d.matrix, 
                                                   len_beta = len_beta)))
    modes <- c(modes, v[which.max(logdensity_glf(par = params, v = v, d.matrix = d.matrix, 
                                                 len_beta = len_beta))])
  }
  pred <- list(mode = modes, density = densities, x = v, formula = formula, newdata = newdata)
  class(pred) <- "predict.ocm"
  return(pred)
}

#' @title Print the output of predict method
#' @description Print method for class \code{predict.ocm}
#' @param x an object of class \code{predict.ocm}
#' @param ... further arguments passed to or from other methods
#' @keywords predict
#' @details The table of predictions from \code{predict.ocm} is printed.
#' @seealso \code{\link{predict.ocm}}, \code{\link{ocm}}
#' @author Maurizio Manuguerra, Gillian Heller
#' @examples 
#' ANZ0001.ocm <- ANZ0001[ANZ0001$cycleno==0 | ANZ0001$cycleno==5,]
#' ANZ0001.ocm$cycleno[ANZ0001.ocm$cycleno==5] <- 1
#' fit.overall <- ocm(overall ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' pred <- predict(fit.overall)
#' print(pred)
#' @export

print.predict.ocm <- function(x, ...)
{
  cat("\nThe data set used by the predict method contains",length(x$mode),"records.\n")
  cat("Call:\n")
  print(update(x$formula, .~.+1))
  cat("\nSummary of modes:\n")
  print(summary(x$mode), ...)
}

#' @title Plot  probability densities  from  output of  predict method
#' @description Plot method for class \code{predict.ocm}
#' @param x An object of class \code{predict.ocm}
#' @param records An integer or a vector of integers. The number of the record/s 
#' in the data set for which the density has to be plotted. If not specified, the 
#' function will  plot all records.
#' @param ... further arguments passed to or from other methods
#' @details The probability densities from \code{predict.ocm}  are plotted.
#' @seealso \code{\link{predict.ocm}}, \code{\link{ocm}}
#' @keywords predict, plot
#' @examples 
#' ANZ0001.ocm <- ANZ0001[ANZ0001$cycleno==0 | ANZ0001$cycleno==5,]
#' ANZ0001.ocm$cycleno[ANZ0001.ocm$cycleno==5] <- 1
#' fit.overall <- ocm(overall ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' pred <- predict(fit.overall)
#' plot(pred)
#' @export
#' @author Maurizio Manuguerra, Gillian Heller

plot.predict.ocm <- function(x, records=NULL, ...)
{
  if (is.null(records)) records=1:nrow(x$density)
  cat("Call:\n")
  print(x$formula)
  cat("The data set used in the predict methos contains ",nrow(x$density)," records.\n")
  for (i in records){
    input <- readline(paste("Press 'enter' to plot the probability density of record ",i,", 'q' to quit: ",sep=''))
    if (input == "q") break()
    plot(x$x, exp(x$density[i,]), ylab="Probability Density", main=paste("Record", i), 
         xlab=paste("mode =", round(x$mode[i],3)), t='l')
    lines(rep(x$mode[i],2), c(0, max(exp(x$density[i,]))), lty=21)
  }
}


#' @title Plot method for Continuous Ordinal Fits
#' 
#' @description Plots the g function as fitted in an \code{ocm} call.
#' @param x an object of class \code{ocm}
#' @param CIs method used for confidence bands for the g function. \code{"no"} = no CIS [default]; \code{"vcov"} = Wald; 
#' \code{"rnd.x.bootstrap"} = random-x bootstrap; \code{"fix.x.bootstrap"} = bootstrap with fixed-x 
#' resampling; \code{"param.bootstrap"} = parametric bootstrap 
#' @param R the number of bootstrap replicates. Ignored if CIs=\code{"no"}
#' @param main  title of the plot. Defauts to ``g function (95\% CIs)''
#' @param xlab  label of the x axis. Defaults to ``Continuous ordinal scale''
#' @param ylab  label of the \code{y} axis. Defaults to an emtpy string
#' @param CIcol  color of the confidence interval bands. Defaults to ``lightblue''
#' @param ... further arguments passed to or from other methods
#' @details The fitted g function of an \code{ocm} object is plotted. 
#' If \code{CIs} is not \code{"no"}, 95\% confidence bands are also plotted.
#' Confidence bands computed with any of the bootstrapping options are 
#' obtained with simple percentiles. 
#' @keywords plot
#' @export
#' @import boot
#' @seealso \code{\link{ocm}}
#' @examples
#' ANZ0001.ocm <- ANZ0001[ANZ0001$cycleno==0 | ANZ0001$cycleno==5,]
#' ANZ0001.ocm$cycleno[ANZ0001.ocm$cycleno==5] <- 1
#' fit.overall  <- ocm(overall  ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' plot(fit.overall, CIs="vcov")
#' \dontrun{
#' plot(fit.overall, CIs="rnd.x.bootstrap", R=100)
#' plot(fit.overall, CIs="fix.x.bootstrap", R=100)
#' plot(fit.overall, CIs="param.bootstrap", R=100)
#' }
#' @author Maurizio Manuguerra, Gillian Heller

plot.ocm <- function(x, CIs = c('no', 'vcov','rnd.x.bootstrap','fix.x.bootstrap','param.bootstrap'), R = 1000, 
                     main="g function (95% CIs)", xlab="Continuous ordinal scale", ylab="", 
                     CIcol = 'lightblue', ...)
{
  #FIXME: this works for glf only: make general?
  #FIXME: with bootstrapping, when a variable is a factor, it could go out of observations for some 
  #level making optim fail? need to use droplevels()
  CIs <- match.arg(CIs)
  R <- as.integer(R)
  len_beta <- x$len_beta
  indices = c(len_beta+1, len_beta+2, len_beta+3)
  params_g <- coef(x)[indices]
  v <- seq(0.01, 0.99, by=0.01)
  gfun <- g_glf(v, params_g)
  xlim <- c(0,1)
  ylim <- c(min(gfun), max(gfun))
  if (CIs=='vcov'){
    #require(MASS)
    vcov_g <- x$vcov[indices, indices]
    #rparams <- mvrnorm(R, params_g, vcov_g, empirical=TRUE)
    rparams <- mvrnormR(R, params_g, vcov_g)
    #FIXME write efficiently
    all_gfuns <- NULL
    for (i in 1:R) all_gfuns <- rbind(all_gfuns, g_glf(v, rparams[i,]))
    ci_low  <- apply(all_gfuns, 2, function(x)quantile(x, 0.025))
    ci_median <- apply(all_gfuns, 2, function(x)quantile(x, 0.5))
    ci_high <- apply(all_gfuns, 2, function(x)quantile(x, 0.975)) 
    ylim <- c(min(ci_low), max(ci_high))
  } else if (CIs=='rnd.x.bootstrap' | CIs=='fix.x.bootstrap'| CIs=='param.bootstrap'){
    #require(boot)
    bs <- boot(x$data, eval(parse(text=CIs)), R, fit = x)
    all_gfuns <- NULL
    for (i in 1:R){
      all_gfuns <- rbind(all_gfuns, g_glf(v, bs$t[i,indices]))
    }
    ci_low  <- apply(all_gfuns, 2, function(x)quantile(x, 0.025))
    #ci_median <- apply(all_gfuns, 2, function(x)quantile(x, 0.5))
    ci_high <- apply(all_gfuns, 2, function(x)quantile(x, 0.975)) 
    ylim <- c(min(ci_low), max(ci_high))
  }
  plot(v, gfun, main=main, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, t='l')
  #CIs
  if (CIs != 'no'){
    #lines(v, ci_low, lty = 2)
    #lines(v, ci_high, lty = 2)
    polygon(c(v, rev(v)),c(ci_low,rev(ci_high)), col = CIcol)
    lines(v,gfun) #to superimpose gfun estimate on shaded area
    #if (CIs=='vcov' | CIs=='rnd.x.bootstrap' | CIs=='fix.x.bootstrap') lines(v, ci_median, lty = 2)
  }
  lines(c(.5,.5), ylim, col='grey')
  lines(xlim, c(0, 0), col='grey')
}

#' @title Anova method for Continuous Ordinal Fits 
#' @description Comparison of continuous ordinal models using likelihood ratio tests. 
#' @param object an object of class \code{ocm}
#' @param ... one or more additional \code{ocm} objects
#' @details Likelihood ratio testing of nested models is performed. 
#' @method anova ocm
#' @keywords anova
#' @export
#' @author Maurizio Manuguerra, Gillian Heller
#'  @seealso \code{\link{ocm}}, \code{\link{print.anova.ocm}}
#' @return The method returns an object of class \code{anova.ocm} and \code{data.frame}, reporting for each model, in hierarchical order:
#' \item{no.par}{number of parameters}
#' \item{AIC}{Akaike information criterion}
#' \item{loglik}{log-likelihood}
#' \item{LR.stat}{likelihood ratio statistic}
#' \item{df}{difference in the degrees of freedom in the models being compared}
#' \item{Pr(>Chisq)}{p-value from the likelihood ratio test}
#' @examples
#' ANZ0001.ocm <- ANZ0001[ANZ0001$cycleno==0 | ANZ0001$cycleno==5,]
#' ANZ0001.ocm$cycleno[ANZ0001.ocm$cycleno==5] <- 1
#' fit.overall  <- ocm(overall  ~ cycleno + bsa + treatment, data=ANZ0001.ocm)
#' anova(fit.overall, update(fit.overall, .~. + age))


anova.ocm <- function(object, ...)
  ### requires that ocm objects have components:
  ###  no.pars: no. parameters used
  ###  call$formula
  ###  link (character)
  ###  gfun (character)
  ###  logLik
  ###
{
  mc <- match.call()
  dots <- list(...)
  ## remove 'test' and 'type' arguments from dots-list:
  not.keep <- which(names(dots) %in% c("test", "type"))
  if(length(not.keep)) {
    message("'test' and 'type' arguments ignored in anova.ocm\n")
    dots <- dots[-not.keep]
  }
  if(length(dots) == 0)
    stop('anova is not implemented for a single "ocm" object')
  mlist <- c(list(object), dots)
  if(!all(sapply(mlist, function(model)
    inherits(model, c("ocm", "ocmm")))))
    stop("only 'ocm' and 'ocmm' objects are allowed")
  nfitted <- sapply(mlist, function(x) length(x$fitted.values))
  if(any(nfitted != nfitted[1L]))
    stop("models were not all fitted to the same dataset")
  no.par <- sapply(mlist, function(x) x$no.pars)
  ## order list with increasing no. par:
  ord <- order(no.par, decreasing=FALSE)
  mlist <- mlist[ord]
  no.par <- no.par[ord]
  no.tests <- length(mlist)
  ## extract formulas, links, gfun:
  forms <- sapply(mlist, function(x) deparse(x$call$formula))
  links <- sapply(mlist, function(x) x$link)
  gfun <- sapply(mlist, function(x) x$gfun)
  models <- data.frame(forms)
  models.names <- c('formula', "link", "gfun")
  models <- cbind(models, data.frame(links, gfun))
  ## extract AIC, logLik, statistics, df, p-values:
  AIC <- sapply(mlist, function(x) -2*x$logLik + 2*x$no.pars)
  logLiks <- sapply(mlist, function(x) x$logLik)
  statistic <- c(NA, 2*diff(sapply(mlist, function(x) x$logLik)))
  df <- c(NA, diff(no.par))
  pval <- c(NA, pchisq(statistic[-1], df[-1], lower.tail=FALSE))
  pval[!is.na(df) & df==0] <- NA
  ## collect results in data.frames:
  tab <- data.frame(no.par, AIC, logLiks, statistic, df, pval)
  tab.names <- c("no.par", "AIC", "logLik", "LR.stat", "df",
                 "Pr(>Chisq)")
  colnames(tab) <- tab.names
  #mnames <- sapply(as.list(mc), deparse)[-1]
  #rownames(tab) <- rownames(models) <- mnames[ord]
  rownames(tab) <- rownames(models) <- paste("Model ",1:length(mlist),":",sep='')
  colnames(models) <- models.names
  attr(tab, "models") <- models
  attr(tab, "heading") <-
    "Likelihood ratio tests of ordinal regression models for continuous scales:\n"
  class(tab) <- c("anova.ocm", "data.frame")
  tab
}


#' @title Print anova.ocm objects
#' 
#' @description Print the results of the comparison of continuous ordinal models in likelihood ratio tests.
#' @param x an object of class \code{anova.ocm}
#' @param digits controls the number of digits to print. Defaults to the maximum of the value 
#' returned by (getOption("digits") - 2) and 3
#' @param signif.stars a logical. Should the significance stars be printed? Defaults to the value 
#' returned by getOption("show.signif.stars")
#' @param ... further arguments passed to or from other methods
#' @keywords summary, anova
#' @seealso \code{\link{ocm}}, \code{\link{anova.ocm}}
#' @return Prints \code{anova.ocm} object
#' @examples
#' ANZ0001.ocm <- ANZ0001[ANZ0001$cycleno==0 | ANZ0001$cycleno==5,]
#' ANZ0001.ocm$cycleno[ANZ0001.ocm$cycleno==5] <- 1
#' fit.overall  <- ocm(overall  ~ cycleno + bsa + treatment, data=ANZ0001.ocm)
#' anova(fit.overall, update(fit.overall, .~. + age))
#' @author Maurizio Manuguerra, Gillian Heller
#' @export

print.anova.ocm <- function(x, digits=max(getOption("digits") - 2, 3), 
                            signif.stars=getOption("show.signif.stars"), ...){
    if (!is.null(heading <- attr(x, "heading")))
      cat(heading, "\n")
    models <- attr(x, "models")
    #row.names(models) <- paste("Model ",1:nrow(models),":",sep='')
    print(models, right=FALSE)
    cat("\n")
    printCoefmat(x, digits=digits, signif.stars=signif.stars,
                 tst.ind=4, cs.ind=NULL, # zap.ind=2, #c(1,5),
                 P.values=TRUE, has.Pvalue=TRUE, na.print="", ...)
    return(invisible(x))
  }



#' @title Extract Log-likelihood for a Continuous Ordinal  Model
#' @description Extracts the log-likelihood for a fitted \code{ocm} object
#' @param object an \code{ocm} object
#' @param ... further arguments to be passed to methods
#' @usage \method{logLik}{ocm}(object, ...)
#' @method logLik ocm
#' @seealso \code{\link{ocm}}
#' @return The log-likelihood of an \code{ocm} object. This is a number with attributes
#' \item{df}{estimated degrees of freedom for the fitted model \code{object}}
#' \item{nobs}{number of observations used in the fitted model \code{object}}
#' \item{class}{class of the returned object: \code{logLik.ocm}}
#' @export
#' @author Maurizio Manuguerra, Gillian Heller
#' @examples
#' ANZ0001.ocm <- ANZ0001[ANZ0001$cycleno==0 | ANZ0001$cycleno==5,]
#' ANZ0001.ocm$cycleno[ANZ0001.ocm$cycleno==5] <- 1
#' fit.overall  <- ocm(overall  ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' logLik(fit.overall)

logLik.ocm <- function(object, ...){
  structure(object$logLik, df = object$df, nobs=object$nobs, class = "logLik.ocm")
}

#' @title Extract AIC from a fitted Continuous Ordinal Model
#' @description Extracts the AIC for a fitted \code{ocm} object
#' @param fit \code{ocm} object
#' @param scale parameter currently not used. For compatibility with general extractAIC method.
#' @param k  ``weight'' of the equivalent degrees of freedom (=: edf) 
#'  in the AIC formula. Defaults to 2
#' @param ... further arguments to be passed to methods
#' @details The generalized AIC is computed:
#' \deqn{-2\ell +k\cdot edf}
#' where \eqn{\ell} is the log likelihood, k=2 gives the AIC, and 
#' k=log(n) gives the BIC.
#' @seealso \code{\link{ocm}}, \code{\link{extractAIC.ocmm}}
#' @return A numeric vector of length 2, with first and second elements giving
#' \item{edf}{the ``equivalent degrees of freedom'' for the fitted model \code{fit}}
#' \item{AIC}{the generalized AIC of \code{ocm} object \code{fit}}
#' @references  Akaike, H (1983). 
#' Information measures and model selection, 
#' \emph{Bulletin of the International Statistical Institute}, 50:277-290.
#' @export
#' @author Maurizio Manuguerra, Gillian Heller
#' @method extractAIC ocm
#' @examples
#' ANZ0001.ocm <- ANZ0001[ANZ0001$cycleno==0 | ANZ0001$cycleno==5,]
#' ANZ0001.ocm$cycleno[ANZ0001.ocm$cycleno==5] <- 1
#' fit.overall  <- ocm(overall  ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' extractAIC(fit.overall)

extractAIC.ocm <- function(fit, scale = 0, k = 2, ...) {
  edf <- fit$df
  c(edf, -2*fit$logLik + k * edf)
}


#' @title Variance-Covariance Matrix for a Fitted Model Object
#' @description Calculates variance-covariance matrix for a fitted \code{ocm} object
#' @param object an \code{ocm} object
#' @param ... further arguments to be passed to methods
#' @details For the generalized logistic g-function, the variance-covariance matrix of model 
#' parameters is 
#' of dimension (\code{len_beta} +3)x(\code{len_beta} +3), where \code{len_beta}  is the number 
#' of beta coefficients in the model.
#' @export
#' @method vcov ocm
#'  @return Variance-covariance matrix of model parameters
#' @seealso \code{\link{ocm}}
#' @author Maurizio Manuguerra, Gillian Heller
#' @examples
#' ANZ0001.ocm <- ANZ0001[ANZ0001$cycleno==0 | ANZ0001$cycleno==5,]
#' ANZ0001.ocm$cycleno[ANZ0001.ocm$cycleno==5] <- 1
#' fit.overall  <- ocm(overall  ~ cycleno + age + bsa + treatment, data=ANZ0001.ocm)
#' vcov(fit.overall)

vcov.ocm <- function(object, ...) {
  object$vcov
}


