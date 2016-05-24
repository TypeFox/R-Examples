
anova.glmrob <- function(object, ...,
                         test = c("Wald", "QD", "QDapprox"))
{
    dotargs <- list(...)
    if (!is.null(names(dotargs))) {
	named <- (names(dotargs) != "")
	if (any(named)) {
	    warning("the following arguments to 'anova.glmrob' are invalid and",
		    "dropped:\n",
		    pasteK(deparse(dotargs[named])))
	    dotargs <- dotargs[!named]
	}
    }
    is.glmrob <- vapply(dotargs, inherits, NA, what="glmrob")
    if(!all(is.glmrob) || !inherits(object, "glmrob"))
	stop("anova.glmrob() only works for 'glmrob' objects")
    test <- match.arg(test)
    if (length(dotargs) > 0)
	anovaGlmrobList(c(list(object), dotargs), test=test)
    else {
	##
	## "'Anova Table' for a single model object
	stop("'Anova Table' for a single model object not yet implemented")
    }
}


anovaGlmrobList <- function (object, test=NULL)
{
  nmodels <- length(object)
  stopifnot(nmodels >= 2)
  responses <- as.character(lapply(object,
				   function(x) deparse(formula(x)[[2]])))
  if (!all(responses == responses[1]))
    stop("Not the same response used in the fitted models")
  nobs <- sapply(object, function(x) length(x$residuals))
  if (any(nobs != nobs[1]))
        stop("models were not all fitted to the same size of dataset")
  methods <- as.character(lapply(object, function(x) x$method))
  if(!all(methods == methods[1]))
    stop("Not the same method used for fitting the models")
  note <- paste("Models fitted by method '", methods[1], "'", sep="")
  tccs <- sapply(object, function(x) length(x$tcc))
  if(!all(tccs == tccs[1]))
    stop("Not the same tuning constant c used in the robust fits")
  ##
  tbl <- matrix(rep(NA, nmodels*4), ncol = 4)
  tbl[1,1] <- nobs[1] - length(coef(object[[1]]))
  for(k in 2:nmodels)
    tbl[k,] <- anovaGlmrobPair(object[[k-1]], object[[k]], test=test)

  ## return
  dimnames(tbl) <- list(1:nmodels,
                        c("pseudoDf", "Test.Stat", "Df", "Pr(>chisq)"))
  title <- switch(test,
                  Wald = "Robust Wald Test Table",
                  QD = "Robust Quasi-Deviance Table",
                  QDapprox =
              "Robust Quasi-Deviance Table Based on a Quadratic Approximation",
                   "")
  variables <- lapply(object, function(x)
                      paste(deparse(formula(x)), collapse = "\n"))
  topnote <- paste("Model ", format(1:nmodels), ": ", variables,
                   sep = "", collapse = "\n")
  structure(as.data.frame(tbl), heading = c(title, "", topnote, note,""),
            class = c("anova", "data.frame"))
}


anovaGlmrobPair <- function(obj1, obj2, test)
{
  if(length(coef(obj1)) < length(coef(obj2))){
      Sign <- 1
      full.mfit <- obj2
      reduced.mfit <- obj1
  }
  else {
      Sign <- -1
      full.mfit <- obj1
      reduced.mfit <- obj2
  }
  X <- model.matrix(full.mfit)
  asgn <- attr(X, "assign")
  tt <- terms(full.mfit)
  tt0 <- terms(reduced.mfit)
  tl <- attr(tt, "term.labels")
  tl0 <- attr(tt0, "term.labels")
  numtl0 <- match(tl0 , tl, nomatch = -1)
  if(attr(tt0, "intercept") == 1) numtl0 <- c(0, numtl0)
  if(any(is.na(match(numtl0, unique(asgn)))))
    stop("Models are not nested!")
  mod0 <- seq(along = asgn)[!is.na(match(asgn, numtl0))]
  if (length(asgn) == length(mod0))
    stop("Models are not strictly nested")
  H0ind <- setdiff(seq(along = asgn), mod0)
  H0coef <- coef(full.mfit)[H0ind]
  df <- length(H0coef)
  pp <- df + length(mod0)

  if(test == "Wald") {
    t.cov <- full.mfit$cov
    t.chisq <- sum(H0coef * solve(t.cov[H0ind, H0ind], H0coef))
    statistic <- c(chisq = t.chisq)
  }
  else if(full.mfit$method=="Mqle" && (test == "QD" || test == "QDapprox")) {
      matM <- full.mfit$matM
      if(test == "QDapprox") {
        ## Difference of robust quasi-deviances
        ## via the asymptotically equivalent quadratic form
        matM11 <- matM[mod0,   mod0, drop=FALSE]
        matM12 <- matM[mod0,  H0ind, drop=FALSE]
        matM22 <- matM[H0ind, H0ind, drop=FALSE]
        matM22.1 <- matM22 - crossprod(matM12, solve(matM11, matM12))
        Dquasi.dev <- nrow(X) * c(H0coef %*% matM22.1 %*% H0coef)
      }
      else {
        quasiDev <- switch(full.mfit$family$family,
			   poisson  = glmrobMqleDiffQuasiDevPois,
			   binomial = glmrobMqleDiffQuasiDevB,
			   Gamma    = glmrobMqleDiffQuasiDevGamma,
                           stop("This family is not implemented"))

        ## note that qdev and qdev0 do depend on an incorrectly specified
        ## lower limits in the integration. But this does't matter in
        ## the following difference, because the difference does not
        ## deepend on it! (Hence I could use the centered nui
        ## (cnui= nui - Enui) in quasiDev as the function to be integrated.
        Dquasi.dev <- quasiDev(mu = full.mfit$fitted.values,
                               mu0 = reduced.mfit$fitted.values,
                               y = full.mfit$y, ni = full.mfit$ni,
                               w.x = full.mfit$w.x, phi=full.mfit$dispersion,
                               tcc = full.mfit$tcc)
      }
      ## Asymptotic distribution: variance and weights of the sum of chi2
      matQ <- full.mfit$matQ
      matM11inv <- solve(matM[mod0,mod0])
      Mplus <- matrix(0, ncol = pp, nrow = pp)
      Mplus[mod0, mod0] <- matM11inv

      d.ev <- Re(eigen(matQ %*% (solve(matM)-Mplus), only.values=TRUE)$values)
      d.ev <- d.ev[1:df] ## just the q (=df) lagest eigenvalues are needed

      if(any(d.ev < 0)) warning("some eigenvalues are negative")

      ## p-value: exact computation for q=1, approximated for q>1 (q=df)
      statistic <- c(quasi.dev = Dquasi.dev/mean(d.ev))

    } else stop("non-implemented test method:", test,
                "for fitting method", full.mfit$method)

  ## return
  c(nrow(X)-pp+df*(Sign<0), Sign*statistic, Sign*df,
    pchisq(as.vector(statistic), df=df, lower.tail = FALSE))
}

