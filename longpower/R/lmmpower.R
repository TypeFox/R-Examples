lmmpower <- function(object, ...) UseMethod("lmmpower")
setGeneric("lmmpower")

lmmpower.default <- function(object=NULL,
   n=NULL,
   parameter = 2,
   pct.change = NULL,
   delta = NULL,
   t = NULL,
   sig.level = 0.05,
   power = NULL, 
   alternative = c("two.sided", "one.sided"),
   beta=NULL,
   beta.CI=NULL,
   delta.CI=NULL,
   sig2.i=NULL,
   sig2.s=NULL,
   sig2.e=NULL,
   cov.s.i=NULL,
   R=NULL,
   method = c("edland", "diggle", "liuliang"),
   ...)
{
	if(sum(!sapply(list(delta, pct.change), is.null))==2) 	
		stop("Only one of delta and pct.change must be specified.")
	if(is.null(delta)&!is.null(beta)&!is.null(pct.change))
	  delta<-pct.change*beta
  if (sum(sapply(list(n, delta, power, sig.level), is.null)) != 1) 
      stop("exactly one of 'n', 'delta', 'power', and 'sig.level' must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
      sig.level | sig.level > 1)) 
      stop("'sig.level' must be numeric in [0, 1]")
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  
	m <- length(t)
	if(is.null(R) & method %in% c("edland", "liuliang")){
	  D <- matrix(c(sig2.i, cov.s.i, cov.s.i, sig2.s), nrow=2)
	  R <- cbind(1,t)%*%D%*%rbind(1,t)  
	  R <- R + diag(sig2.e, m, m)
	}

	if(method == "liuliang"){
	  u <- list(u1 = t, u2 = rep(0,m))
	  v <- list(v1 = cbind(1,1,rep(0,m)),
	         v2 = cbind(1,0,t))
	  if(!is.null(n)) N <- n*2 else N <- NULL
	}

	results <- switch(method,
	  edland = edland.linear.power(n=n, delta=delta, t=t, 
      sig2.s=sig2.s, sig2.e=sig2.e, 
      sig.level=sig.level,
      power=power,
      alternative=alternative,...),
	  diggle = diggle.linear.power(n=n, delta=delta, t=t, R=R, 
	    sig.level=sig.level,
	    power=power,
	    alternative=alternative,...),
	  liuliang = liu.liang.linear.power(N=N, delta=delta, u=u, v=v, R=R,
	    sig.level=sig.level,
	    power=power,
	    alternative=alternative,...))

	if(is.null(delta.CI)&!is.null(beta.CI)) results$delta.CI <- (results$delta/beta)*beta.CI
	if(!is.null(beta)) results$beta <- beta
	if(!is.null(beta.CI)) results$beta.CI <- beta.CI
		
	if(!is.null(results$delta.CI)){
		n.upper <- switch(method,
		  edland = edland.linear.power(n=NULL, results$delta.CI[1], t=t, sig2.s, sig2.e, 
  	    sig.level=sig.level,
  	    power=power,
  	    alternative=alternative,...)$n,
      diggle = diggle.linear.power(n=NULL, results$delta.CI[1], t=t, R=R, 
		    sig.level=sig.level,
		    power=power,...)$n,
		  liuliang = liu.liang.linear.power(N=NULL, results$delta.CI[1], u=u, v=v, R=R, 
		    sig.level=sig.level,
		    power=power,...)$N/2)
		n.lower <- switch(method,
		  edland = edland.linear.power(n=NULL, results$delta.CI[2], t, sig2.s, sig2.e, 
        sig.level=sig.level,
        power=power,
        alternative=alternative,...)$n,
      diggle = diggle.linear.power(n=NULL, results$delta.CI[2], t=t, R=R, 
		    sig.level=sig.level,
		    power=power,...)$n,
		  liuliang = liu.liang.linear.power(N=NULL, results$delta.CI[2], u=u, v=v, R=R, 
		    sig.level=sig.level,
		    power=power,...)$N/2)
		n.CI <- c(n.lower, n.upper)
		if(n.CI[1]>n.CI[2]) n.CI <- n.CI[2:1]
		results$n.CI <- n.CI 
	}

    if(is.character(parameter)){
      names(results)[names(results) == "beta"] <- parameter
      names(results)[names(results) == "beta.CI"] <- paste(parameter, "CI")
  }
	results <- results[unlist(lapply(results, function(x) !is.null(x)))]
	structure(results, class = "power.longtest")
}

lmmpower.lme <- function(object,
   n = NULL,
   parameter = 2,
   pct.change = NULL,
   delta = NULL,
   t = NULL,
   sig.level = 0.05,
   power = NULL, 
   alternative = c("two.sided", "one.sided"),
   beta=NULL,
   beta.CI=NULL,
   delta.CI=NULL,
   sig2.i=NULL,
   sig2.s=NULL,
   sig2.e=NULL,
   cov.s.i=NULL,
   method = c("edland", "diggle", "liuliang"),
   ...)
{
	alternative <- match.arg(alternative)
  method <- match.arg(method)
  
	if(is.numeric(parameter)) parameter <- rownames(summary(object)$tTable)[parameter]
	
	tab <- nlme::getVarCov(object)

	if(nrow(tab)>2) stop("Too many random effects. Function is 
	  	equipped to handle at most a random intercept and slope.")

	if(is.null(beta))	
	  beta = summary(object)$tTable[parameter,'Value']
	if(is.null(beta.CI))	
	  beta.CI = rep(summary(object)$tTable[parameter,'Value'],2) + 
	    c(-1,1)*qnorm(0.025)*summary(object)$tTable[parameter,'Std.Error']
	if(beta.CI[1] > beta.CI[2]) beta.CI <- beta.CI[2:1]
	# var of random intercept
	if(is.null(sig2.i))
	  sig2.i = tab["(Intercept)", "(Intercept)"]
	# var of random slope
	if(is.null(sig2.s))
	  sig2.s = ifelse(nrow(tab)==1, 0, 
	         ifelse(nrow(tab)==2, tab[2, 2], NA))
	
	# residual var
	if(is.null(sig2.e))
	  sig2.e = object$sigma^2
	# covariance of slope and intercep
	if(is.null(cov.s.i))
	  cov.s.i = ifelse(nrow(tab)==1, 0, 
	          ifelse(nrow(tab)==2, tab[1, 2], NA))

	lmmpower(object=NULL,
		n = n,
		parameter = parameter,
		pct.change = pct.change,
		delta = delta,
		t = t,
		sig.level = sig.level,
		power = power, 
		alternative = alternative,
		beta = beta,
		beta.CI = beta.CI,
		delta.CI = delta.CI,
		sig2.i = sig2.i,
		sig2.s = sig2.s,
		sig2.e = sig2.e,
		cov.s.i = cov.s.i, 
		method = method, ...)
}

lmmpower.gee <- function(object,
   n = NULL,
   parameter = 2,
   pct.change = NULL,
   delta = NULL,
   t = NULL,
   sig.level = 0.05,
   power = NULL, 
   alternative = c("two.sided", "one.sided"),
   beta=NULL,
   beta.CI=NULL,
   delta.CI=NULL,
   method = c("diggle", "liuliang"),
   ...)
{
	alternative <- match.arg(alternative)
  method <- match.arg(method)
  
	if(is.numeric(parameter)) parameter <- rownames(summary(object)$coefficients)[parameter]
	
	if(is.null(beta))	
	  beta = summary(object)$coefficients[parameter,'Estimate']
	if(is.null(beta.CI))	
	  beta.CI = rep(summary(object)$coefficients[parameter,'Estimate'],2) + 
	    c(-1,1)*qnorm(0.025)*summary(object)$coefficients[parameter,'Robust S.E.']
	if(beta.CI[1] > beta.CI[2]) beta.CI <- beta.CI[2:1]
	
	R <- summary(object)$working.correlation * summary(object)$scale

	lmmpower(object=NULL,
		n = n,
		parameter = parameter,
		pct.change = pct.change,
		sig.level = sig.level,
		power = power, 
		alternative = alternative,
		sig2.e=1,
		sig2.s=0,
		beta=beta,
		beta.CI=beta.CI,
		delta.CI=delta.CI,
		R=R,
		t=t, 
		method=method, ...)
}

setMethod("lmmpower", signature(object = "merMod"),
  function(object, 
   n = NULL, 
   parameter = 2,
   pct.change = NULL,
   delta = NULL,
   t = NULL,
   sig.level = 0.05,
   power = NULL, 
   alternative = c("two.sided", "one.sided"),
   beta=NULL,
   beta.CI=NULL,
   delta.CI=NULL,
   sig2.i=NULL,
   sig2.s=NULL,
   sig2.e=NULL,
   cov.s.i=NULL,
   method = c("edland", "diggle", "liuliang"),
   ...)
{
  if (!(sum(sapply(list(n, delta, power, sig.level), is.null)) == 1 |
        sum(sapply(list(n, pct.change, power, sig.level), is.null)) == 1)) 
      stop("exactly one of 'n', 'delta' (or 'pct.change'), 'power', and 'sig.level' must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
      sig.level | sig.level > 1)) 
      stop("'sig.level' must be numeric in [0, 1]")
	alternative <- match.arg(alternative)
  method <- match.arg(method)
  
	if(is.numeric(parameter)) parameter <- rownames(coef(summary(object)))[parameter]
	
	tab <- lme4::VarCorr(object)
	if(length(tab)>1) stop("Too many grouping levels. Function is 
	  	equipped to handle at most one grouping level.")
	tab <- tab[[1]]

	if(nrow(tab)>3) stop("Too many random effects. Function is 
	  	equipped to handle at most a random intercept and slope.")
  m <- length(t)
	if(is.null(beta))	
	  beta = fixef(object)[parameter][[1]]
	if(is.null(beta.CI))	
	  beta.CI = rep(coef(summary(object))[parameter,"Estimate"],2) + 
	    c(-1,1)*qnorm(0.025)*coef(summary(object))[parameter,"Std. Error"]
	if(beta.CI[1] > beta.CI[2]) beta.CI <- beta.CI[2:1]

	# var of random intercept
	if(is.null(sig2.i))
	  sig2.i = tab[1, 1]
	# var of random slope
	if(is.null(sig2.s))
	  sig2.s = ifelse(nrow(tab)==1, 0, 
	         ifelse(nrow(tab)==2, tab[2 ,2], NA))
	# residual var
	if(is.null(sig2.e))
	  sig2.e = getME(object, "sigma")^2
	# covariance of slope and intercept
	if(is.null(cov.s.i))
    cov.s.i = ifelse(nrow(tab)==1, 0, 
            ifelse(nrow(tab)==2, tab[2, 1], NA))

	lmmpower(n=n, object=NULL,
		parameter = parameter,
		pct.change = pct.change,
		delta = delta,
		t = t,
		sig.level = sig.level,
		power = power, 
		alternative = alternative,
		beta=beta,
		beta.CI=beta.CI,
		delta.CI=delta.CI,
		sig2.i=sig2.i,
		sig2.s=sig2.s,
		sig2.e=sig2.e,
		cov.s.i=cov.s.i, 
		method = method, ...)
})
