"tukey.1df" <- function(aov.obj, data, error.term="Within")
{
  multistratum <- inherits(aov.obj, what="aovlist")
  if (multistratum)
  {  res.df <- aov.obj[[error.term]]$df.residual
     aov.proj <- proj(aov.obj)
     res <- aov.proj[[error.term]][,"Residuals"]
     fit <- fitted.errors(aov.obj, error.term = error.term)
     rno1 <- attr(attr(aov.obj, "terms"), "response")+1
     rname <- as.character(attr(attr(aov.obj, "terms"), "variables"))[rno1]
  }
  else
  {	res.df <- aov.obj$df.residual
    res <- residuals(aov.obj)
    fit <- fitted(aov.obj)
    rno1 <- attr(attr(aov.obj, "terms"), "response")+1
    rname <- as.character(attr(attr(aov.obj, "terms"), "variables"))[rno1]
  }
  data.temp <- data
  data.temp[,rname] <- fit*fit
  if (multistratum)
  {	Tukey.aov <- aov(attr(aov.obj, "terms"), data.temp)
    Tukey.proj <- proj(Tukey.aov)
    res2 <- Tukey.proj[[error.term]][,"Residuals"]
    res2
  }
  else
  {	Tukey.aov <- aov(aov.obj$terms, data.temp)
    res2 <- residuals(Tukey.aov)
  }
  SS.1df.num <- sum(res * res2)^2
  SS.1df.den <- sum(res2^2)
  SS.1df <- SS.1df.num/SS.1df.den
  if (SS.1df.den < 1e-10) cat("** Warning - there appears to be extremely little non-linear variation so that","\n",
                              "  the values for Tukey.SS are unstable and the results below may be unreliable.","\n",
                              "  Only use if at least two non-interacting factors above the same Residual","\n",
                              "  in the analysis.","\n")
  SS.res <- sum(res^2)
  SS.res.1 <- SS.res - SS.1df
  F.1df <- ((res.df - 1) * SS.1df)/SS.res.1
  p.1df <- 1 - pf(F.1df, 1, res.df-1)
  return(list(Tukey.SS = SS.1df, Tukey.F = F.1df, Tukey.p = p.1df, Devn.SS = SS.res.1))
}

"residuals.aovlist" <- function(object, error.term=NULL, ...)
{
	if (!inherits(object, what="aovlist"))
	  stop("Not an aovlist object")
  res <- 0
	aov.proj <- proj(object)
	if (is.null(error.term))
	  res <- aov.proj[[length(object)]][,"Residuals"]
	else
		res <- aov.proj[[error.term]][,"Residuals"]
	res
}

"fitted.aovlist" <- function(object, error.term=NULL, ...)
{
	if (!inherits(object, what="aovlist"))
	  stop("Not an aovlist object")
	aov.proj <- proj(object)
	if (is.null(error.term))
		no.strata <- length(object)
	else
		no.strata <- which(names(aov.proj) == error.term)
	fit <- aov.proj[["(Intercept)"]][,1]
	for (i in 2:no.strata)
	{	nterms <- ncol(aov.proj[[i]])
		if (dimnames(aov.proj[[i]])[[2]][nterms] == "Residuals") nterms <- nterms-1
		if (nterms > 0)
			if (nterms == 1)
				fit <- fit + aov.proj[[i]][,1]
			else
				fit <- fit + rowSums(aov.proj[[i]][,1:nterms])
	}
  fit
}

resid.errors <- function(...) residuals.aovlist(...)
fitted.errors <- function(object, error.term=NULL, ...) 
                               fitted.aovlist(object=object, error.term=error.term, ...)

#"resid.errors" = function(object, ...) UseMethod("residuals")
#fitted.errors = function(x) UseMethod("fitted.errors")
#setGeneric("fitted.errors.aovlist")
