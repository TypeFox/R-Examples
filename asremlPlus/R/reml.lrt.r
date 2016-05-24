reml.lrt.asreml <- function(full.asreml.obj, 
                            reduced.asreml.obj, 
                            positive.zero=FALSE)
{
#asreml codes:
#  B - fixed at a boundary (!GP)
#  F - fixed by user
#  ? - liable to change from P to B    
#  P - positive definite
#  C - Constrained by user (!VCC)      
#  U - unbounded
#  S - Singular Information matrix

  #Check have asreml objectst
  if ((is.null(full.asreml.obj) | is.null(reduced.asreml.obj)) |
        class(full.asreml.obj) != "asreml" | 
        class(reduced.asreml.obj) != "asreml")
    stop("Must supply two objects of class 'asreml'")

  #Check that fixed and sparse models are the same
  fixed.labels <- lapply(list(full.asreml.obj,reduced.asreml.obj), 
                         function(x) {attr(terms(x$fixed.formula), "term.labels")})
  sparse.labels <- lapply(list(full.asreml.obj,reduced.asreml.obj), 
                          function(x) {attr(terms(x$sparse), "term.labels")})
  mu <- sapply(list(full.asreml.obj,reduced.asreml.obj), 
               function(x) {attr(terms(x$fixed.formula), "intercept")})
  if (!all(mu == mu[1])) 
    stop("Fixed models must be identical with repsect to the intercept")
  if (!all(diff(sapply(fixed.labels, function(x) length(x))) == 0)) 
    stop("Fixed models differ in length")
  if (all(sapply(fixed.labels, function(x) length(x)) > 0))
  { for (i in 2:length(fixed.labels)) 
    { if (!all(is.element(fixed.labels[[1]], fixed.labels[[i]]))) 
        stop("Fixed models differ")
    }
  }
  if (!all(diff(sapply(sparse.labels, function(x) length(x))) == 0)) 
    stop("sparse models differ in length")
  if (all(sapply(sparse.labels, function(x) length(x)) > 0)) 
  { for (i in 2:length(sparse.labels)) 
    { if (!all(is.element(sparse.labels[[1]], sparse.labels[[i]]))) 
        stop("sparse models differ")
    }
  }
  
  #Perform the test
  summ.full <- summary(full.asreml.obj)
	summ.reduce <- summary(reduced.asreml.obj)
	REMLLRT <- 2*(summ.full$loglik-summ.reduce$loglik)
	vc.full <- summ.full$varcomp
	vc.reduce <- summary(reduced.asreml.obj)$varcomp
	DF <- nrow(vc.full) - nrow(vc.reduce)
	if ("Fixed" %in% levels(vc.full$constraint))
		DF <- DF - table(vc.full$constraint)["Fixed"]
	if ("Constrained" %in% levels(vc.full$constraint))
		DF <- DF - table(vc.full$constraint)["Constrained"]
	if ("Singular" %in% levels(vc.full$constraint))
	  DF <- DF - table(vc.full$constraint)["Singular"]
	if  ("Fixed" %in% levels(vc.reduce$constraint))
		DF <- DF + table(vc.reduce$constraint)["Fixed"]
	if  ("Constrained" %in% levels(vc.reduce$constraint))
		DF <- DF + table(vc.reduce$constraint)["Constrained"]
	if ("Singular" %in% levels(vc.reduce$constraint))
	  DF <- DF + table(vc.reduce$constraint)["Singular"]
	DF <- unname(DF)
	if (DF == 0)
	  warning("DF is zero indicating no difference between models in the number of parameters")
   else
     if (DF <= 0)
       warning("Negative degrees of freedom indicating the second model is more complex")
   if (positive.zero)
   { if (REMLLRT < 1e-10)
	 	 p <- 1
	  else if (DF == 1)
#this option is used when testing a positively-contrained variance component is zero
#it adjusts for testing on the boundary of the contraint by computing 0.5(ch1_0 + chi_2)
#comes from Self & Liang (1987) Case 5
        p <- 0.5*(1-pchisq(REMLLRT, DF))
  	 else
#following is for DF positively-contrained components equal to zero 
#computes sum_i=0 to DF mix_i * chi_i where mix_i = choose(DF,  i)*2^(_DF)
#comes from Self & Liang (1987) Case 9 with s = 0
     {  df <- seq(DF)
        p <- 1-2^(-DF)-sum(choose(DF, df)*2^(-DF)*pchisq(REMLLRT, df))
	  }
   }
   else
	     p <- 1-pchisq(REMLLRT, DF)
   data.frame(REMLLRT, DF, p)
}


info.crit.asreml <- function(asreml.obj)
{
	summ <- summary(asreml.obj)
	vc <- summ$varcomp
	DF <- nrow(vc)
	if ("Fixed" %in% levels(vc$constraint))
		DF <- DF - table(vc$constraint)["Fixed"]
	if ("Constrained" %in% levels(vc$constraint))
		DF <- DF - table(vc$constraint)["Constrained"]
	if ("Singular" %in% levels(vc$constraint))
	  DF <- DF - table(vc$constraint)["Singular"]
	names(DF) <- ""
	logREML <- summ$loglik
#calculate AIC and BIC
	AIC <- -2 * logREML + 2 * DF
	BIC <- -2 * logREML + DF * log(summ$nedf)
	data.frame(DF, AIC, BIC, logREML)
}

#reml.lrt <- function(full.asreml.obj, reduced.asreml.obj, 
#                     positive.zero=FALSE)
#{ 
#  test <- reml.lrt.asreml(full.asreml.obj = full.asreml.obj, 
#                          reduced.asreml.obj = reduced.asreml.obj, 
#                          positive.zero=positive.zero)
#  return(test)
#}

#info.crit <- function(asreml.obj)
#{ info <- info.crit.asreml(asreml.obj = asreml.obj)
#  return(info)
#}

reml.lrt <- function(...)
{ .Deprecated(new = "reml.lrt.asreml", package = "asremlPlus")
  invisible()
}

info.crit <- function(...)
{ .Deprecated(new = "info.crit.asreml", package = "asremlPlus")
  invisible()
}
