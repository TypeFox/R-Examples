dispersiontest <- function(object, trafo = NULL, alternative = c("greater", "two.sided", "less"))
{
  if(!inherits(object, "glm") || family(object)$family != "poisson")
    stop("only Poisson GLMs can be tested")
  alternative <- match.arg(alternative)
  otrafo <- trafo
  if(is.numeric(otrafo)) trafo <- function(x) x^otrafo
  
  y <- if(is.null(object$y)) model.response(model.frame(object)) else object$y
  yhat <- fitted(object)
  aux <- ((y - yhat)^2 - y)/yhat
    
  if(is.null(trafo)) {
    STAT <- sqrt(length(aux)) * mean(aux)/sd(aux)
    NVAL <- c(dispersion = 1)
    EST <- c(dispersion = mean(aux) + 1)    
  } else {
    auxreg <- lm(aux ~ 0 + I(trafo(yhat)/yhat))
    STAT <- as.vector(summary(auxreg)$coef[1,3])
    NVAL <- c(alpha = 0)
    EST <- c(alpha = as.vector(coef(auxreg)[1]))
  }
  
  rval <- list(statistic = c(z = STAT),
               p.value = switch(alternative,
                           "greater" = pnorm(STAT, lower.tail = FALSE),
		           "two.sided" = pnorm(abs(STAT), lower.tail = FALSE)*2,
		           "less" = pnorm(STAT)),
	       estimate = EST,
	       null.value = NVAL,
	       alternative = alternative,
	       method = switch(alternative,
	                  "greater" = "Overdispersion test",
			  "two.sided" = "Dispersion test",
			  "less" = "Underdispersion test"),
	       data.name = deparse(substitute(object)))
  class(rval) <- "htest"
  return(rval)
}


## NB. score tests a la DCluster now implemented in countreg
##
## TODO:
## LRT for Poi vs NB2.
## fix DCluster::test.nb.pois() and pscl::odTest()
## proposed interface:
##   poistest(object, object2 = NULL)
## where either a "negbin" and a "glm" object have to be
## supplied or only one of them, then update via either
##   cl <- object$call
##   cl[[1]] <- as.name("glm.nb")
##   cl$link <- object$family$link
##   cl$family <- NULL
## or
##   cl <- object$call
##   cl[[1]] <- as.name("glm")
##   cl$family <- call("poisson")
##   cl$family$link <- object$family$link
##   cl$link <- NULL
##   cl$init.theta <- NULL
## and evaluate the call "cl" appropriately.
