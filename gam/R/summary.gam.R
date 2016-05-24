"summary.gam" <-
   function (object, dispersion = NULL, ...) 
{
  object.lm=object
  class(object.lm)="lm"
  paod=anova(object.lm,...)
  attr(paod,"heading")="Anova for Parametric Effects"

  save.na.action <- object$na.action
  object$na.action <- NULL
  fun <- function(assign, coeff) sum(!is.na(coeff[assign]))
  wt <- object$weights
  coef <- object$coef
  dresid <- residuals(object, "deviance")
  resid <- object$residuals
  n <- length(resid)
  s <- object$s
  nl.chisq <- object$nl.chisq
  assg <- object$assign
  if (is.null(assg)) 
    assg <- attributes(object$terms)$assign
  df <- rep(1, length(assg))
  df[is.na(object$coef)] <- 0
  df <- tapply(df, assg, sum)
  dfnames <- attr(object$terms, "term.labels")
  if (attr(object$terms, "intercept")) 
    dfnames <- c("(Intercept)", dfnames)
  names(df) <- dfnames
  df <- unlist(df)
  nldf <- object$nl.df
  n <- length(object$residuals)
  if (is.null(rdf <- object$df.resid)) {
    rdf <- n - sum(df)
    if (!is.null(nldf)) 
      rdf <- rdf - sum(nldf)
  }
  if (!is.null(wt)) {
    wt <- wt^0.5
    resid <- resid * wt
    excl <- wt == 0
    if (any(excl)) {
      warning(paste(sum(excl), "rows with zero weights not counted"))
      resid <- resid[!excl]
      dresid <- dresid[!excl]
      if (is.null(object$df.residual)) 
        rdf <- rdf - sum(excl)
    }
  }
  if (rdf > 0) 
    phihat <- sum(resid^2)/rdf
  else {
    phihat <- Inf
    warning("Residual degrees of freedom are negative or zero.  This occurs when the sum of the parametric and nonparametric degrees of freedom exceeds the number of observations.  The model is probably too complex for the amount of data available.")
  }
  famname <- object$family[["family"]]
  if (is.null(famname)) 
    famname <- "gaussian"
  chiorf <- TRUE
  if (!is.null(dispersion) && dispersion == 0) 
    dispersion <- phihat
  if (is.null(dispersion)) 
    dispersion <- switch(famname, poisson = 1, binomial = 1, 
                         {
                           chiorf <- FALSE
                           phihat
                         })
  names(dispersion) <- famname
  if (length(nldf)) {
    aod <- as.matrix(round(df, 1))
    dimnames(aod) <- list(names(df), "Df")
    if (!is.null(nl.chisq)) {
      aod <- cbind(aod, NA, NA, NA)
      nl.chisq <- nl.chisq/dispersion
      snames <- names(nldf)
      aod[snames, 2] <- round(nldf, 1)
      aod[snames, 3] <- if (!chiorf) 
        nl.chisq/nldf
      else nl.chisq
      aod[snames, 4] <- if (chiorf) 
        1 - pchisq(nl.chisq, nldf)
      else if (rdf > 0) 
        1 - pf(nl.chisq/nldf, nldf, rdf)
      else NA
      rnames <- c("Df", "Npar Df", "Npar Chisq", "P(Chi)")
      if (!chiorf) 
        rnames[3:4] <- c("Npar F", "Pr(F)")
      dimnames(aod) <- list(names(df), rnames)
      heading <- "Anova for Nonparametric Effects"
    }
    else heading <- "DF for Nonparametric Terms"
    aod <- as.anova(data.frame(aod[,-1], check.names = FALSE), 
                    heading)
  }
  else aod <- NULL
  structure(list(call = object$call, terms = object$terms, 
                 anova = aod, parametric.anova=paod, dispersion = dispersion, df = c(sum(df) + 
                                                         sum(nldf), rdf), deviance.resid = dresid, deviance = deviance(object), 
                 null.deviance = object$null.deviance, aic = object$aic, 
                 iter = object$iter, na.action = save.na.action), class = "summary.gam")
}
