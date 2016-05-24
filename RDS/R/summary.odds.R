#' Summarizing Generalized Linear Model Fits with Odds Ratios for Survey Data
#' 
#' \code{RDS::summary.svyglm.RDS} is a version of \code{summary.svyglm} that 
#' reports odds-ratios in place of coefficients in the summary table. 
#' This only applies for the \code{binomial} family. Otherwise it is identical to
#' \code{summary.svyglm}.
#' The default in \code{summary.svyglm} is to display the log-odds-ratios 
#' and this displays the exponetiated from
#' and a 95% confidence interval in place of the standard errors and \code{z ratio} columns. The
#' p-values are still displayed.
#' 
#' \code{svyglm} fits a generalised linear model to data from a complex survey design, with
#' inverse-probability weighting and design-based standard errors.
#' 
#' There is no \code{anova} method for \code{svyglm} as the models are not
#' fitted by maximum likelihood.
#' 
#' See the manual page on \code{svyglm} for detail of that function.
#' @param object an object of class \code{"svyglm"}, usually, a result of a call
#' to \code{\link[survey]{svyglm}}.
#' @param correlation logical; if \code{TRUE}, the correlation matrix of the
#' estimated parameters is returned and printed.
#' @param df.resid Optional denominator degrees of freedom for Wald tests.
#' @param odds logical; Should the coefficients be reported as odds (rather than log-odds)?
#' @param \dots further arguments passed to or from other methods.
#' @return \code{RDS::summary.svyglm} returns an object of class \code{"summary.svyglm.RDS"},
#' a list with components
#' 
#' \item{call}{the component from \code{object}.} \item{family}{the component
#' from \code{object}.} \item{deviance}{the component from \code{object}.}
#' \item{contrasts}{the component from \code{object}.} \item{df.residual}{the
#' component from \code{object}.} \item{null.deviance}{the component from
#' \code{object}.} \item{df.null}{the component from \code{object}.}
#' \item{deviance.resid}{the deviance residuals: see
#' \code{\link[survey]{residuals.svyglm}}.} \item{coefficients}{the matrix of
#' coefficients, standard errors, z-values and p-values.  Aliased coefficients
#' are omitted.} \item{aliased}{named logical vector showing if the original
#' coefficients are aliased.} \item{dispersion}{either the supplied argument or
#' the inferred/estimated dispersion if the latter is \code{NULL}.} \item{df}{a
#' 3-vector of the rank of the model and the number of residual degrees of
#' freedom, plus number of coefficients (including aliased ones).}
#' \item{cov.unscaled}{the unscaled (\code{dispersion = 1}) estimated
#' covariance matrix of the estimated coefficients.} \item{cov.scaled}{ditto,
#' scaled by \code{dispersion}.} \item{correlation}{(only if \code{correlation}
#' is true.)  The estimated correlations of the estimated coefficients.}
#' \item{symbolic.cor}{(only if \code{correlation} is true.)  The value of the
#' argument \code{symbolic.cor}.}
#' \item{odds}{Are the coefficients reported as odds (rather than log-odds)?}
#' @seealso \code{\link[survey]{svyglm}}, \code{\link{summary}}.
#' @keywords models regression
#' @examples
#' 
#' ## For examples see example(svyglm)
#' 
#' @export
#' @method summary svyglm.RDS
summary.svyglm.RDS<-function (object, correlation = FALSE, df.resid=NULL, odds=TRUE, ...) 
{
    Qr <- object$qr
    est.disp <- TRUE
    if (is.null(df.resid))
      df.r <- object$df.residual
    else
      df.r<-df.resid
    
    dispersion<-survey::svyvar(stats::resid(object,"pearson"), object$survey.design,
                       na.rm=TRUE)
    
    coef.p <- stats::coef(object)
    covmat<-stats::vcov(object)
    dimnames(covmat) <- list(names(coef.p), names(coef.p))
    var.cf <- diag(covmat)
    s.err <- sqrt(var.cf)
    tvalue <- coef.p/s.err
    if(object$family$family == "binomial" & odds){
     dn <- c("Estimate", "Lower")
     qvalue <- stats::qt(0.975,df=df.r)
     if (!est.disp) {
        pvalue <- 2 * pnorm(-abs(tvalue))
        coef.table <- cbind(exp(coef.p), 
          exp(coef.p-qvalue*s.err), 
	  exp(coef.p+qvalue*s.err), pvalue)
        dimnames(coef.table) <- list(names(coef.p),
                                         c(dn, "Upper","Pr(>|z|)"))
     }
     else if (df.r > 0) {
        pvalue <- 2 * stats::pt(-abs(tvalue), df.r)
        coef.table <- cbind(exp(coef.p), 
              exp(coef.p-qvalue*s.err), 
	      exp(coef.p+qvalue*s.err), pvalue)
        dimnames(coef.table) <- list(names(coef.p),
                                     c(dn, "Upper","Pr(>|t|)"))
     }
     else {
        coef.table <- cbind(exp(coef.p), Inf)
        dimnames(coef.table) <- list(names(coef.p), dn)
     }
    }else{
     dn <- c("Estimate", "Std. Error")
     if (!est.disp) {
        pvalue <- 2 * pnorm(-abs(tvalue))
        coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
        dimnames(coef.table) <- list(names(coef.p), c(dn, "z value", 
            "Pr(>|z|)"))
     }
     else if (df.r > 0) {
        pvalue <- 2 * stats::pt(-abs(tvalue), df.r)
        coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
        dimnames(coef.table) <- list(names(coef.p), c(dn, "t value", 
            "Pr(>|t|)"))
     }
     else {
        coef.table <- cbind(coef.p, Inf)
        dimnames(coef.table) <- list(names(coef.p), dn)
     }
    }
    ans <- c(object[c("call", "terms", "family", "deviance", 
        "aic", "contrasts", "df.residual", "null.deviance", "df.null", 
        "iter")], list(deviance.resid = stats::residuals(object, type = "deviance"), 
        aic = object$aic, coefficients = coef.table, dispersion = dispersion, 
        df = c(object$rank, df.r,NCOL(Qr$qr)), cov.unscaled = covmat, 
        cov.scaled = covmat))
    if (correlation) {
        dd <- sqrt(diag(covmat))
        ans$correlation <- covmat/outer(dd, dd)
    }
    ans$aliased<-is.na(stats::coef(object,na.rm=FALSE))
    ans$survey.design<-list(call=object$survey.design$call)
    ans$odds<-odds
    class(ans) <- c("summary.svyglm.RDS","summary.svyglm","summary.glm")
    return(ans)
}

#' Summarizing Generalized Linear Model Fits with Odds Ratios
#' 
#' \code{print.summary.svyglm.RDS} is a version of \code{print.summary.svyglm} that 
#' reports odds-ratios in place of coefficients in the summary table. 
#' This only applies for the \code{binomial} family. Otherwise it is identical to
#' \code{print.summary.svyglm}.
#' The default in\cr
#' \code{print.summary.svyglm} is to display the log-odds-ratios 
#' and this displays the exponetiated from
#' and a 95% confidence interval in place of the standard errors and \code{z ratio} columns. The
#' p-values are still displayed.
#' 
#' @param x an object of class \code{"summary.svyglm.RDS"}, usually, a result of a
#' call to \code{RDS::summary.svyglm}.
#' @param digits the number of significant digits to use when printing.
#' @param symbolic.cor logical. If \code{TRUE}, print the correlations in a
#' symbolic form (see \code{\link{symnum}}) rather than as numbers.
#' @param signif.stars logical. If \code{TRUE}, \sQuote{significance stars} are
#' printed for each coefficient.
#' @param \dots further arguments passed to or from other methods.
#' @seealso \code{\link[survey]{svyglm}}, \code{\link[survey]{summary.svyglm}}.
#' @keywords models regression
#' @examples
#' 
#' ## For examples see example(svyglm)
#' 
#' @export
#' @method print summary.svyglm.RDS
print.summary.svyglm.RDS<-function (x, digits = max(3, getOption("digits") - 3),
                                symbolic.cor = x$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...) 
{
  ##if (!exists("printCoefmat")) stats::printCoefmat<-print.coefmat

  cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "") 

    cat("Survey design:\n")
    print(x$survey.design$call)
   
        if (!is.null(df <- x$df) && (nsingular <- df[3] - df[1])) 
         if(x$odds & x$family$family == "binomial"){
            cat("\nOdds Ratios: (", nsingular, " not defined because of singularities)\n", 
                sep = "")
	 }else{
            cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", 
                sep = "")
	 }
        else {
        if(x$odds & x$family$family == "binomial"){
          cat("\nOdds Ratios\n")
	}else{
          cat("\nCoefficients\n")
	}
	}
        coefs <- x$coefficients
        if (!is.null(aliased <- is.na(x$coefficients[,1])) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn, 
                colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
        }
#       if(x$family$family == "binomial"){
#  coefs[,1:3] <- exp(coefs[,1:3])
#}
	stats::printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
    
    cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ", 
        format(x$dispersion), ")\n\n",  "Number of Fisher Scoring iterations: ", 
        x$iter, "\n", sep = "")
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            if (is.logical(symbolic.cor) && symbolic.cor) {
                print(stats::symnum(correl, abbr.colnames = NULL))
            }
            else {
                correl <- format(round(correl, 2), nsmall = 2, 
                  digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop = FALSE], quote = FALSE)
            }
        }
    }
    cat("\n Warning: The model fit currently only adjusts for the first-order RDS weights, and not other complexities of the sampling.\n  The numerical summaries should be interpreted with caution.\n")
    cat("\n")
    invisible(x)
}
