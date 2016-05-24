#################################################################################
##
##   R package parma
##   Alexios Ghalanos Copyright (C) 2012-2013 (<=Aug)
##   Alexios Ghalanos and Bernhard Pfaff Copyright (C) 2013- (>Aug)
##   This file is part of the R package parma.
##
##   The R package parma is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package parma is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################
# Shannon Entropy
sentropy = function(w){
	if(any(w)<0){
		idx = which(w<0)
		if(all(abs(w[idx])<1e-4)) w[idx] = abs(w[idx]) else stop("\nShannon Entropy only valid for positive weights")
	}
	return( -sum(w * log(w)) )
}

# Cross Entropy Approximation of Golan, Judge and Miller (1996)
# w = reference weights
# q = portfolio weights
centropy = function(w, q){
	if(any(q)<0){
		idx = which(q<0)
		if(all(abs(q[idx])<1e-4)) q[idx] = abs(q[idx]) else stop("\nCross Entropy weights q must be positive")
	}
	return( sum( (1/q) * (w - q)^2 ) )
}


vfun.mad = function(Data){
	return( abs(scale(Data, scale = F)) )
}

vfun.lpm = function(Data, threshold, moment){
	if(threshold == 999){
		Data = scale(Data, scale = FALSE)
		threshold =  0
	}
	return( pmax(Data - threshold, 0)^moment )
}

vfun.var = function(Data){
	return( scale(Data, scale = FALSE)^2 )
}

vfun.cvar = function(Data, VaR = NULL, alpha){
	if(is.null(VaR)){
		n = dim(Data)[1]
		Rp = Data
		sorted = sort(Rp)
		n.alpha = floor(n * alpha)
		VaR = sorted[n.alpha]
	}
	d = VaR - Data
	f = -VaR + func.max.smooth(d)
	return( f )
}


fun.mad = function(w, Data, benchmark = NULL){
	d = Data %*% w
	if(is.null(benchmark)) d = scale(d, scale=FALSE) else d = d - benchmark
	return( mean(abs(d)) )
}

fun.lpm = function(w, Data, threshold, moment){
	if(threshold == 999){
		Data = scale(Data, scale = FALSE)
		threshold =  0
	}
	return( (mean(pmax(Data %*% w - threshold, 0)^moment))^(1/moment) )
}

fun.upm = function(w, Data, threshold, moment){
	if(threshold == 999){
		Data = scale(Data, scale = FALSE)
		threshold =  0
	}
	return( (mean(pmax(threshold - Data %*% w, 0)^moment))^(1/moment) )
}

fun.minmax = function(w, Data, benchmark = NULL){
	if(is.null(benchmark)) benchmark = 0
	return( min(Data %*% w - benchmark) )
}

fun.var = function(w, Data, benchmark = NULL){
	d = Data %*% w
	if(is.null(benchmark)) d = scale(d, scale=FALSE) else d = d - benchmark	
	return( mean(d^2) )
}

fun.qvar = function(w, S){
	return( w %*% S %*% w)
}

fun.cvar = function(w, Data, VaR = NULL, alpha = 0.05, benchmark = NULL){
	if(is.null(benchmark)) benchmark = 0
	if(is.null(VaR)){
		n = dim(Data)[1]
		Rp = (Data %*% w - benchmark)
		sorted = sort(Rp)
		n.alpha = floor(n * alpha)
		VaR = sorted[n.alpha]
	}
	d = VaR - (Data %*% w - benchmark)
	f = -VaR + mean( func.max.smooth(d) )/alpha
	return( f )
}

# DrawDown and RunUp
fun.cumret = function(w, Data, benchmark = NULL)
{
	if(is.null(benchmark)) benchmark = 0
	f = cumsum( as.numeric( Data %*% w - benchmark) )
	return( f )
}

fun.dd = function(w, Data, lower = TRUE, benchmark = NULL)
{
	if(is.null(benchmark)) benchmark = 0
	# Absolute Drawdown
	x = fun.cumret(w, Data, benchmark)
	n = length(x)
	if( lower ){
		f = sapply( 1:n, FUN = function(i) max(x[1:i]) - x[i] )
	} else{
		f = sapply( 1:n, FUN = function(i) x[i] - min(x[1:i]) )
	}
	return( f )
}

fun.cdar = function(w, Data, DaR = NULL, alpha = 0.05, lower = TRUE, benchmark = NULL)
{	
	if(is.null(benchmark)) benchmark = 0
	port = Data %*% w - benchmark
	if( !lower ) port = - port
	dd = fun.dd(w, Data, lower, benchmark)
	n = length( dd )
	probability = rep( 1/n, n )
	xdd = sort.int(dd, index.return = TRUE)
	sidx  = xdd$ix
	# key step: reverse so that we work with 0:alpha (left tail)
	xdd = rev(xdd$x)
	cp = cumsum(probability[sidx[1:n]])
	intalpha = max( max( which( cp <= alpha ) ), 1)
	if(is.null(DaR)) DaR = as.numeric( xdd[intalpha] )
	
	intalpha2 = max( max( which( cp < alpha ) ), 1 )
	intalpha3 = max( min( which( cp >= alpha ) ), 1)
	
	CDaRplus = sum( probability[sidx[1:intalpha2]]/sum(probability[sidx[1:intalpha2]]) * xdd[1:intalpha2] )
	lambda = ( sum(probability[1:intalpha3]) - alpha )/( alpha )
	f = as.numeric( lambda * DaR + (1 - lambda) * CDaRplus )
	
	return( f )
}


fun.rachevratio = function(w, Data, alphadn, alphaup){
	dncvar = fun.cvar(w,  Data, VaR = NULL, alpha = alphadn, benchmark = 0)
	upcvar = fun.cvar(w, -Data, VaR = NULL, alpha = alphaup, benchmark = 0)
	return( dncvar/upcvar )
}


riskfun = function(weights, Data, risk = c("mad", "ev", "minimax", "cvar", "cdar", "lpm"), 
		benchmark = NULL, alpha = 0.05, moment = 1, threshold = 0, VaR = NULL, DaR = NULL)
{
	xrisk = match.arg(tolower(risk[1]), c("mad", "ev", "minimax", "cvar", "cdar", "lpm"))	
	ans = switch(tolower(xrisk),
			mad = fun.mad(weights, Data, benchmark),
			ev = fun.var(weights, Data, benchmark),
			minimax = fun.minmax(weights, Data, benchmark),
			cvar = fun.cvar(weights, Data, VaR = VaR, alpha = alpha, benchmark),
			lpm = fun.lpm(weights, Data, moment = moment, threshold = threshold),
			cdar = fun.cdar(weights, Data, DaR = DaR, alpha = alpha, lower = TRUE, benchmark))
	return(abs(ans))
}

fun.risk = function(weights, Data, options, risk = c("mad", "ev", "minimax", "cvar", "cdar", "lpm"), 
		benchmark = NULL, VaR = NULL, DaR = NULL)
{
	xrisk = match.arg(tolower(risk[1]), c("mad", "ev", "minimax", "cvar", "cdar", "lpm"))	
	ans = switch(xrisk,
			mad = fun.mad(weights, Data, benchmark),
			ev = fun.var(weights, Data, benchmark),
			minimax = fun.minmax(weights, Data, benchmark),
			cvar = fun.cvar(weights, Data, VaR = VaR, alpha = options$alpha, benchmark),
			lpm = fun.lpm(weights, Data, moment = options$moment, threshold = options$threshold),
			cdar = fun.cdar(weights, Data, DaR = DaR, alpha = options$alpha, lower = TRUE, benchmark))
	return(abs(ans))
}


vfun.risk = function(scenario, options, risk){
	ans = switch(tolower(risk),
			mad = vfun.mad(scenario),
			ev = vfun.var(scenario),
			cvar = vfun.cvar(scenario, VaR = NULL, alpha = options$alpha),
			lpm = vfun.lpm(scenario, moment = options$moment, threshold = options$threshold))
	return(ans)
}


simweights = function(m, LB = 0, UB = 1, budget=1, forecast=rep(0,m), constrain.positive = TRUE, rseed = NULL){
	if(is.null(rseed)) rseed = as.integer(Sys.time())
	bin = rbinom(m, 1, pnorm(rnorm(1)))
	i=0
	wx = rep(0, m)
	while(i==0){
		nb = sum(bin)
		w = rtruncnorm(nb, a=LB, b=UB, mean = (UB-LB)/2, sd = (UB-LB)/3)
		w = budget * (w/sum(w))
		wx[as.logical(bin)] = w
		if(constrain.positive){
			if(sum(wx*forecast)>=0) i=1
		} else{
			i=1
		}
	}
	return(wx)
}

is.even = function(x){ as.logical( (round(x)+1)%%2 ) }

.logtransform = function(x, LB, UB){ LB + (UB - LB)/(1+exp(-x))}

# matrix square root via svd
.sqrtm = function(x)
{
	tmp = svd(x)
	sqrtx = tmp$u%*%sqrt(diag(tmp$d))%*%t(tmp$u)
	return(sqrtx)
}
