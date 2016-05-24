logconcave = function(times, covariates, aug = TRUE){
  if(missing(covariates))
    return(logconcave.univariate(times) )
  recenter = FALSE
  if(class(times) == "numeric"){
    n <- length(times)
    augVal = 0
    if(aug == TRUE)
      augVal = 1/(2*n)
    
    x <- unique(sort(times))
    
    max_ind = which(times == max(times))
    if(length(max_ind) > 1)
      stop('Problem: two uncensored values tied for last. Need to decide what to do about this. See Anderson-Bergman 2013 appendix E')
    
    covariates <- as.matrix(covariates)
 #   center_Covar <- covariates[max_ind, ]
 #   covariates <- t(t(covariates) - center_Covar)
    
    if(max(x) == Inf | min(x) == -Inf)
      stop("Error: non finite values supplied")
    b = rep(0, length(x))
    output = .Call('LC_CoxPH', times, times, x, b, as.integer( c(1, length(x) ) ), FALSE, augVal, augVal, as.matrix(covariates) )
    names(output) <- c("x", "phi", "llk", "beta", "full_Hess")
   	l_beta = length(output$beta)
   	hess_dim = length(output$full_Hess[,1])
    output$beta_Covariance = -solve(output$full_Hess[-hess_dim, -hess_dim])[1:l_beta, 1:l_beta]
    output[['dens']] <- exp(output$phi)
    if(recenter == FALSE)
      output[['covOffset']] <- rep(0, length(output$beta))
    class(output) <- "LCPHObject"
    return(output)
  }
  
#  if(is.matrix(times))
   int.Data = as.matrix(times)
   
   if(ncol(int.Data) != 2)
   	stop("times should either be a vector of length n or be of dimensions n x 2")
    
  l = int.Data[,1]
  u = int.Data[,2]
  if(max(l) <= min(u) )
  {
    x.use = c(max(l), min(u))
    density = c(1/(min(u) - max(l) ),1/(min(u) - max(l) ))
    lk.final = 0
    it = 0
    x = x.use
    l.dens = log(density)
    stop("Estimator degenerate; universal overlap in data")
    #		output <- list(x.use, density, lk.final, it, x, l.dens)
    #	cat("Note: Universal overlap leads to poorly performing estimator!\n")
    #	names(output) <- c("x", "dens", "lk", "it", "x.all", "l.dens")
    #	return(output)
  }	
  
  
  output = start.Inf(l, u)
  x = output[[1]]
  b = output[[2]]
  
  max_time = max(u)
  max_ind = which(l == max_time)
  if(length(max_ind) == 0)
    max_ind = 1
  if(length(max_ind) > 1)
    stop('Problem: two uncensored values tied for last. Need to decide what to do about this. See Anderson-Bergman 2013 appendix E')
  
  covariates <- as.matrix(covariates)
  
#  center_Covar <- covariates[max_ind, ]
#  covariates <- t(t(covariates) - center_Covar)
  n <- length(l)
  augVal = 0
  if(aug == TRUE)
    augVal = 1/(2*n)
  minX = if(x[1] == -Inf) 2 else 1
  maxX = if(x[length(x)] == Inf) length(x) - 1 else length(x) 
  actInds = c(minX, which(b == max(b) ) , maxX )
  inds <- make.inds(x, l,u)
  output = .Call('LC_CoxPH', l, u, x, b, as.integer(actInds), TRUE, augVal, augVal, as.matrix(covariates) )
  names(output) <- c("x", "phi", "llk", "beta", "full_Hess")
  l_beta = length(output$beta)
  hess_dim = length(output$full_Hess[,1])
  output$beta_Covariance = -solve(output$full_Hess[-hess_dim, -hess_dim])[1:l_beta, 1:l_beta]
  output$beta_Covariance = try(-solve(output$full_Hess)[1:l_beta, 1:l_beta])
  output[['dens']] <- exp(output$phi)
  if(recenter == FALSE)
    output[['covOffset']] <- rep(0, length(output$beta))
  class(output) <- "LCPHObject"
  return(output)
}

simPH_Censored <- function(n = 100, b1 = 0.5, b2 = -0.5, shape = 2){
  rawQ <- runif(n)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  nu <- exp(x1 * b1 + x2 * b2)
  trueTime <- rgamma(n, rate = nu, shape = shape)
  cTime <- rgamma(n, rate = 1, shape = shape)
  leftCensored <- cTime > trueTime
  left <- rep(0,n)
  right <- rep(Inf, n)
  left[!leftCensored] <- cTime[!leftCensored]
  right[leftCensored] <- cTime[leftCensored]
  output <- list(times = cbind(left, right), x = cbind(x1, x2))
  return(output)
}

sim_Censored <- function(n = 100){
	l = rep(0, n)
	u = rep(1, n)
	time = rbeta(n, 2, 2)
	cTime = runif(n)
	isLeftCens <- time < cTime
	u[isLeftCens] <- cTime[isLeftCens]
	l[!isLeftCens] <- cTime[!isLeftCens]
	return(cbind(l, u))
}


logconcave.univariate = function(data){
	if(class(data) == "numeric"){
		x = sort(data)
		repVec <- as.numeric(table(x))
		x <- unique(x)
		inds = 1:length(x)
		if(max(x) == Inf | min(x) == -Inf)
			stop("Error: non finite values supplied")
		b = rep(0, length(x))
		output = .Call('uniVarLCCens', inds, inds, x, b, repVec, as.integer( c(1, length(x) ) ), FALSE )
		names(output) <- c("x", "beta", "llk", 'iterations')
		output[['dens']] <- exp(output$beta)
		class(output) <- "LCObject"
		return(output)
	}

	int.Data = as.matrix(data)
	if(ncol(data) != 2)
		stop("times should either be a matrix of length n or be of dimensions n x 2")
	l = int.Data[,1]
	u = int.Data[,2]
	

	
	if(max(l) <= min(u) )
		{
		x.use = c(max(l), min(u))
		density = c(1/(min(u) - max(l) ),1/(min(u) - max(l) ))
		lk.final = 0
		it = 0
		x = x.use
		l.dens = log(density)
		output <- list(x.use, density, lk.final, it, x, l.dens)
		cat("Note: Universal overlap leads to poorly performing estimator!\n")
		names(output) <- c("x", "dens", "lk", "it", "x.all", "l.dens")
	#	return(output)
		}	

	
	output = start.Inf(l, u)
	x = output[[1]]
	b = output[[2]]
	repData <- rep.vec(int.Data)
	l <- repData[,1]
	u <- repData[,2]
	repVec = repData[,3]
	minX = if(x[1] == -Inf) 2 else 1
	maxX = if(x[length(x)] == Inf) length(x) - 1 else length(x) 
	actInds = c(minX, which(b == max(b) ) , maxX )
	actInds = unique(actInds)
	inds <- make.inds(x, l,u)
	output <- .Call('uniVarLCCens', inds[,'l.inds'], inds[,'u.inds'], x, b, repVec, as.integer(actInds) , TRUE)
	names(output) <- c("x", "beta", "llk", 'iterations')
	output[['dens']] <- exp(output$beta)
	class(output) <- "LCObject"
	return(output)
}


plotLC <- function(fit, funtype = 'surv', covars = NA, ...){
  if(class(fit) == 'LCObject')
    plot_LCObject(fit, funtype = funtype, ...)
  if(class(fit) == 'LCPHObject')
    plot_LCPHObject(fit, funtype = funtype, covars = covars, ...)
}

linesLC <- function(fit, funtype = 'surv', covars = NA, ...){
  if(class(fit) == 'LCObject')
    lines_LCObject(fit, funtype = funtype, covars = covars, ...)
  if(class(fit) == 'LCPHObject')
    lines_LCPHObject(fit, funtype = funtype, covars = covars, ...)
}

plot_LCObject <- function(x, funtype = "surv", ...){
  fit <- x
  x <- fit$x
	x <- x[x != -Inf & x!= Inf]
	minx = min(x)
	maxx = max(x)
	grid = 0:200/200 * (maxx - minx) + minx
	if(funtype == 'surv'){
		y <- 1 - pLC(grid, fit)
#		if(missing(xlab))
			xlab = 't'
#		if(missing(ylab))
			ylab = 'S(t)'
		}
	else if(funtype == 'pdf'){
		y <- dLC(grid, fit)
#		if(missing(xlab))
			xlab = 'x'
#		if(missing(ylab))	
			ylab = 'f(x)'
		}
	else if(funtype == 'cdf'){
		y <- pLC(grid, fit)
#		if(missing(xlab))
			xlab = 'x'
#		if(missing(ylab))
			ylab = 'F(x)'
		}
	else
		stop('funtype not recongized. Valid choices are "surv", "pdf" and "cdf"')
	plot(grid, y, xlab = xlab, ylab = ylab, type = 'l', ... )
}

lines_LCObject <- function(x, y, ...){
  if(missing(y))
    y <- 'surv'
  fit <- x
  funtype <- y
  x <- fit$x
	x <- x[x != -Inf & x!= Inf]
	minx = min(x)
	maxx = max(x)
	grid = 0:200/200 * (maxx - minx) + minx
	if(funtype == 'surv'){
		y <- 1 - pLC(grid, fit)
#		if(missing(xlab))
			xlab = 't'
#		if(missing(ylab))
			ylab = 'S(t)'
		}
	else if(funtype == 'pdf'){
		y <- dLC(grid, fit)
#		if(missing(xlab))
			xlab = 'x'
#		if(missing(ylab))	
			ylab = 'f(x)'
		}
	else if(funtype == 'cdf'){
		y <- pLC(grid, fit)
#		if(missing(xlab))
			xlab = 'x'
#		if(missing(ylab))
			ylab = 'F(x)'
		}
	else
		stop('funtype not recongized. Valid choices are "surv", "pdf" and "cdf"')
	lines(grid, y, xlab = xlab, ylab = ylab, ... )
}

plot_LCPHObject <- function(x, covars, funtype = "surv", ...){
  if(any(is.na(covars)))
    stop('valid covars argument required for CoxPH model')  
  fit <- x
  x <- fit$x
  x <- x[x != -Inf & x!= Inf]
  minx = min(x)
  maxx = max(x)
  grid = 0:200/200 * (maxx - minx) + minx
  if(funtype == 'surv'){
    y <- 1 - pLC(grid, fit, covars)
    #		if(missing(xlab))
    xlab = 't'
    #		if(missing(ylab))
    ylab = 'S(t)'
  }
  else if(funtype == 'pdf'){
    y <- dLC(grid, fit, covars)
    #		if(missing(xlab))
    xlab = 'x'
    #		if(missing(ylab))	
    ylab = 'f(x)'
  }
  else if(funtype == 'cdf'){
    y <- pLC(grid, fit, covars)
    #		if(missing(xlab))
    xlab = 'x'
    #		if(missing(ylab))
    ylab = 'F(x)'
  }
  else
    stop('funtype not recongized. Valid choices are "surv", "pdf" and "cdf"')
  plot(grid, y, xlab = xlab, ylab = ylab, type = 'l', ... )
}
lines_LCPHObject <- function(x, funtype = "surv", covars, ...){
  if(any(is.na(covars)))
    stop('valid covars argument required for CoxPH model')
  fit <- x
  x <- fit$x
  x <- x[x != -Inf & x!= Inf]
  minx = min(x)
  maxx = max(x)
  grid = 0:200/200 * (maxx - minx) + minx
  if(funtype == 'surv'){
    y <- 1 - pLC(grid, fit, covars)
    #		if(missing(xlab))
    xlab = 't'
    #		if(missing(ylab))
    ylab = 'S(t)'
  }
  else if(funtype == 'pdf'){
    y <- dLC(grid, fit, covars)
    #		if(missing(xlab))
    xlab = 'x'
    #		if(missing(ylab))	
    ylab = 'f(x)'
  }
  else if(funtype == 'cdf'){
    y <- pLC(grid, fit, covars)
    #		if(missing(xlab))
    xlab = 'x'
    #		if(missing(ylab))
    ylab = 'F(x)'
  }
  else
    stop('funtype not recongized. Valid choices are "surv", "pdf" and "cdf"')
  lines(grid, y, xlab = xlab, ylab = ylab, ... )
}

makeNonCS_RepVec <- function(l, u){
	ord = order(l)
	sort_L = l[ord]
	sort_R = u[ord]
	out_l = NA
	out_r = NA
	out_rep = NA
	count = 0
	n = length(sort_L)
	for(i in 1:n){
		l_val = sort_L[i]
		if(!is.na(l_val)){
			count = count + 1
			r_val = sort_R[i]
			areEqual = which(sort_L[i:n] == l_val & sort_R[i:n] == r_val)
			out_l[count] <- l_val
			out_r[count] <- r_val
			out_rep[count] <- length(areEqual)
			sort_L[i + areEqual -1 ] <- NA
		}
	}
	return(list(out_l, out_r, out_rep))
}


rep.vec <- function(cs.data)
	{
	drop <- min(cs.data[,1]) == cs.data[,1] & max(cs.data[,2]) == cs.data[,2]
	cs.data <- cs.data[!drop,]
	l <- cs.data[,1]
	u <- cs.data[,2]
	n <- length(l)
	min.l <- min(l)
	max.u <- max(u)
	left <- u[l == min.l]
	right <- l[u == max.u]	
	left <- sort(left)
	right <- sort(right)
	if(length(left) + length(right) != n)
		{
		rep.vec.output <- makeNonCS_RepVec(l, u)
		return(cbind(rep.vec.output[[1]], rep.vec.output[[2]], rep.vec.output[[3]], row.names = NULL) ) 
		}
	rep.lft <- table(left)
	val.lft <- unique(left)
	rep.rgt <- table(right)
	val.rgt <- unique(right)
	l <- c( rep(min.l, length(val.lft)), val.rgt )
	u <- c( val.lft, rep(max.u, length(val.rgt) ) )
	rep.vec <- c(rep.lft, rep.rgt) 
	names(rep.vec) <- NULL
		return(cbind(l, u, rep.vec,  row.names = NULL) )
	}		
	
make.inds <- function(x, l, u)
	{
	u.inds <- NA
	l.inds <- NA
	for(i in 1:length(l) ) 
		{
		u.inds[i] <- which(x == u[i])
		l.inds[i] <- which(x == l[i])
		}
	return(cbind(u.inds, l.inds) ) 
	}
	
	
start.Inf <- function(l, u)
	{
	x.pos <- sort( unique( c(l, u) ) )
	
	k = length(x.pos)
	
	x.mid = (x.pos[-1] + x.pos[-k])/2
	if(x.pos[1] == -Inf)
		x.mid[1] = x.pos[2] - 1
	if(x.pos[k] == Inf)
		x.mid[k-1] = x.pos[k-1] + 1
	x <- sort( c(x.pos, x.mid) )
	k <- length(x)
	
	mid.ind = round(k/2)
	x.mid = x[mid.ind]
	
	
	betas = -abs(x - x.mid)/sd(x[abs(x) < Inf])
	
	return(list(x, betas) )
	}
dLC = function(x, fit, covars)
	{
	if(class(fit) != "LCObject" & class(fit) != "LCPHObject")
		{
		cat("Error: LCObject must be from logconcave function!\n")
		return(NULL)
		}
	if(!is.numeric(x))
		{
		cat("Error: x must be numeric!\n")
		return(NULL)
		}
	dom = fit$x
	phi = log(fit$dens)
	x.dens = rep(0, length(x) ) 
	for(i in 1:length(x) ) 
		x.dens[i] = .Call('d_lc', as.numeric(x[i]), as.numeric(dom), as.numeric(phi) )
#		{
#			{
#			l.ind = max(which(dom <= x[i]) )
#			r.ind = min(which(dom >= x[i]) )
#			if(l.ind == r.ind)
#				x.dens[i] = exp(phi[l.ind])
#			else if(phi[l.ind] > -Inf & phi[r.ind] > -Inf)
#				{
#				est.phi = phi[l.ind] + (phi[r.ind] - phi[l.ind])/(dom[r.ind] - dom[l.ind]) * (x[i] - dom[l.ind])
#				x.dens[i] = exp(est.phi)
#				}
#			else if(phi[l.ind] == -Inf)
#				{
#				slp = (phi[2] - phi[3])/(dom[3] - dom[2])
#				est.phi = phi[2] + slp * (dom[2] - x[i])
#				x.dens[i] = exp(est.phi)
#				}
#			else if(phi[r.ind] == -Inf)
#				{
#				slp = (phi[l.ind] - phi[l.ind-1])/(x[l.ind] - x[l.ind-1])
#				est.phi = phi[l.ind] + slp *(x[i] - dom[l.ind])
#				x.dens[i] = exp(est.phi)
#				}
#			}
#		}
  if(class(fit) == "LCPHObject"){
#    covars <- as.matrix(covars)
    if(length(dim(covars)) < 2){
      numCovars = length(covars)
      covars <- matrix(covars, nrow = 1)
    }
    else
      numCovars <- dim(covars)[2]
    if(numCovars != length(fit$beta) & length(fit$beta) != 1)
      stop( paste('number of covariates provided not equal to number of regression coefficients') )
    s.ests <- 1 - pLC(x, fit, covars * 0)
    adjCovars <- t(t(covars) - fit$covOffset)
    nu <- exp(adjCovars %*% fit$beta)
    x.dens <- as.numeric(nu * x.dens * s.ests^(nu - 1) ) 
  }
	return(x.dens)
	}	
	
pLC = function(x, fit, covars)
	{
	if(class(fit) != "LCObject" & class(fit) != 'LCPHObject')
		{
		cat("Error: LCObject must be from logconcave function!\n")
		return(NULL)
		}
	if(!is.numeric(x))
		{
		cat("Error: x must be numeric!\n")
		return(NULL)
		}
	dom = fit$x
	phi = log(fit$dens)
	x.p = rep(0, length(x) ) 
#	x.p[x >= max(dom)] <- 1
	for(i in 1:length(x) ) 
		x.p[i] = .Call('p_lc', as.numeric(x[i]), as.numeric(dom), as.numeric(phi) )
#		{
#		if(x[i] > min(dom) & x[i] < max(dom) & x[i] >-Inf & x[i] < Inf ) 
#			{
#			l.ind = max(which(dom <= x[i]) )
#			r.ind = min(which(dom >= x[i]) ) 
#			if(l.ind == r.ind)
#				{
#				p.est = 0
#				for(j in 1:(l.ind-1) )
#					p.est = p.est + exact.int(phi[j], phi[j+1], dom[j], dom[j+1])
#				x.p[i] = p.est
#				}
#			if(phi[l.ind] > -Inf & phi[r.ind] > -Inf)
#				{
#				est.phi = phi[l.ind] + (phi[r.ind] - phi[l.ind])/(dom[r.ind] - dom[l.ind]) * (x[i] - dom[l.ind])
#				p.est = 0
#				if(l.ind > 1)
#					{
#					for(j in 1:(l.ind-1) )
#						p.est = p.est + exact.int(phi[j], phi[j+1], dom[j], dom[j+1])
#					}
#				p.est = p.est + exact.int(phi[l.ind], est.phi, dom[l.ind], x[i])
#				x.p[i] = p.est
#				}
#			if(dom[l.ind] == -Inf)
#				{
#        cat('oops...there\'s an error in pLC. Email cianders@uci.edu\n')
#				}
#			if(dom[r.ind] == Inf)
#				{
#				slp = (phi[l.ind] - phi[l.ind-1])/(dom[l.ind] - dom[l.ind-1])
#				est.phi = phi[l.ind] + slp *(x[i] - dom[l.ind])
#				for(j in 1:(l.ind-1) )
#					p.est = p.est + exact.int(phi[j], phi[j+1], dom[j], dom[j+1])
#				p.est = p.est + exact.int(phi[l.ind], est.phi, dom[l.ind], x[i])
#				x.p[i] = p.est
#				}
#			}
#		}
	if(class(fit) == "LCPHObject"){
    if(length(dim(covars)) < 2){
      numCovars <- length(covars)[1]
      covars <- matrix(covars, nrow = 1)
    }
    else
      numCovars <- dim(covars)[2]
	  if(numCovars != length(fit$beta) & length(fit$beta) != 1)
	    stop( paste('number of covariates provided not equal to number of regression coeffecients') )
	          adjCovars <- t(t(covars) - fit$covOffset)
	          nu <- exp(adjCovars %*% fit$beta)
            x.p = as.numeric( (1 - (1 - x.p)^nu) )
	}	
	return(x.p)
	}		

qLC.opt.fxn = function(q, p, fit)
	return( (pLC(q, fit) - p) ^ 2)

qLCCov.opt.fxn = function(q, p, fit, covars)
  return( (pLC(q, fit, covars) - p) ^ 2)


qLC <- function(p, fit, covars)
	{
	q = rep(NA, length(p))
	dom = fit$x[fit$x > -Inf & fit$x < Inf]
	interval = c(min(dom), max(dom) ) 
	if(class(fit) == 'LCObject'){
	  p.interval = pLC(interval, fit)
	  for(i in 1:length(p) )
  		{
  		if(p[i] > p.interval[1] & p[i] < p.interval[2])
  	 	q[i] <- optimize(qLC.opt.fxn, interval = interval, p = p[i], fit = fit)$minimum
  	 	}
	}
  if(class(fit) == 'LCPHObject'){
    if (length(dim(covars)) < 2) {
      numCovars <- length(covars)[1]
      covars <- matrix(covars, nrow = 1)
    }
    p.interval = pLC(interval, fit, covars)
    for(i in 1:length(p) )
    {
      if(p[i] > p.interval[1] & p[i] < p.interval[2])
        q[i] <- optimize(qLCCov.opt.fxn, interval = interval, p = p[i], fit = fit, covars = covars[i,])$minimum
    }
  }
  return(q)
	}
	
exact.int <- function(b.l, b.u, x.l, x.u)
	{
	if(x.l == x.u)
		return(0)
	if(b.l == b.u)
		return(exp(b.l) * (x.u - x.l) ) 
	if(x.l == -Inf)
		{
		slp = (b.l - b.u) / (x.u - x.l)
		return(-exp(b.u) / slp)
		}
	if(x.u == Inf)
		{
		slp = (b.u - b.l) / (x.u - x.l)
		return(-exp(b.l) / slp)
		}
	output = (exp(b.u) - exp(b.l) ) *  (x.u - x.l) /(b.u - b.l) 
	return(output)
	}		
	
	

icpower = function(data, alpha = 2){
	if(class(data) == "numeric"){
		x = sort(data)
		repVec <- as.numeric(table(x))
		x <- unique(x)
		inds = 1:length(x)
		if(max(x) == Inf | min(x) == -Inf)
			stop("Error: non finite values supplied")
		b = rep(1, length(x))
		output = .Call('uniVarICCens', inds, inds, x, b, repVec, as.integer( c(1, length(x) ) ), FALSE, as.numeric(alpha) )
		names(output) <- c("x", "beta", "llk")
		output[['dens']] <- (output$beta)^(-alpha)
		class(output) <- "ICPObject"
		return(output)
	}

if(is.matrix(data))
	int.Data = data
else
	stop("Data type not recognized")	
l = int.Data[,1]
u = int.Data[,2]
	

	
	if(max(l) <= min(u) )
		{
		x.use = c(max(l), min(u))
		density = c(1/(min(u) - max(l) ),1/(min(u) - max(l) ))
		lk.final = 0
		it = 0
		x = x.use
		l.dens = log(density)
		output <- list(x.use, density, lk.final, it, x, l.dens)
		cat("Note: Universal overlap leads to poorly performing estimator!\n")
		names(output) <- c("x", "dens", "lk", "it", "x.all", "l.dens")
	#	return(output)
		}	

	
	
	#NEED TO DECIDE ON WHAT TO DO ABOUT START VALUE FOR INVERSE CONVEX POWER...
	output = start.Inf(l, u)
	x = output[[1]]
	b = output[[2]]
	b = -b + 1
	repData <- rep.vec(int.Data)
	l <- repData[,1]
	u <- repData[,2]
	repVec = repData[,3]
	minX = if(x[1] == -Inf) 2 else 1
	maxX = if(x[length(x)] == Inf) length(x) - 1 else length(x) 
	actInds = c(minX, which(b == min(b) ) , maxX )
	actInds= unique(actInds)
	inds <- make.inds(x, l,u)
	output <- .Call('uniVarICCens', inds[,'l.inds'], inds[,'u.inds'], x, b, repVec, as.integer(actInds) , TRUE, as.numeric(alpha) )
	names(output) <- c("x", "beta", "llk")
	output[['dens']] <- (output$beta)^(-alpha)
	class(output) <- "ICPObject"
	return(output)
}

	
	