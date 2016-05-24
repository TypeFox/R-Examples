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

################################################################
# General Notes:
# wm represents the number of assets
# fm represents the problem size
# widx represents the indices of the weights
# midx is the index of the multiplier for the optimal formulation
# vidx is the index of the VaR for the CVaR problem
################################################################


################################################################
# General Functions, Constants and smooth approximations
################################################################

parma.Small = 1e-20
parma.Big = 1e20
func.max.smooth = function(x)
{
	return( 0.5 * (sqrt(x*x + parma.Small ) + x) )
}

grad.max.smooth = function(x)
{
	return( 0.5 + 0.5*x/sqrt(x*x+parma.Small) )
}

func.abs.smooth = function(x)
{
	return( sqrt( ( x + parma.Small ) * ( x + parma.Small ) ) )
}

grad.abs.smooth = function(x)
{
	return( 0.5*(2*x+2*parma.Small)/sqrt((x+parma.Small)*(x+parma.Small)) )
}
func.abs.smooth2 = function(x)
{
	return( (2 * x/pi) * atan(parma.Big * x) )
}

grad.abs.smooth2 = function(x)
{
	return( (2 * 1/pi) * atan(parma.Big*x) + (2 * x/pi) * 1/(1+parma.Big*x*x) )
}
##################################################################
# Variance (EV)
##################################################################

func.minvar = function(w, optvars, uservars){
	N = optvars$N
	widx = optvars$widx
	f = sum( ( optvars$Data %*% as.numeric(w[widx]) - optvars$benchmark)^2)/N
	return( f )
}

grad.minvar = function(w, optvars, uservars){
	Data = optvars$Data
	wm = optvars$wm
	fm = optvars$fm
	widx = optvars$widx
	N = optvars$N
	g = vector(mode = "numeric", length = wm)
	port = as.numeric( Data %*% as.numeric(w[widx]))
	xport = port - optvars$benchmark
	for(i in 1:wm){
		g[i] =  sum( 2*(Data[,i]*xport) )/N
	}
	fg = rep(0, fm)
	fg[widx] = g
	return( fg )
}


func.optvar = function(w, optvars, uservars){
	N = optvars$N
	widx = optvars$widx
	midx = optvars$midx
	f =  sum( ( (optvars$Data %*% as.numeric(w[widx])) - w[midx]*optvars$benchmark)^2 )/N
	return( f )
}

grad.optvar = function(w, optvars, uservars)
{
	wm = optvars$wm
	fm = optvars$fm
	N = optvars$N
	widx = optvars$widx
	midx = optvars$midx
	Data = optvars$Data
	port = Data %*% w[widx]
	g = vector(mode = "numeric", length = wm)
	xport = port - w[midx]*optvars$benchmark
	for(i in 1:wm){
		g[i] =  sum( 2*(Data[,i]*xport)/N )
	}
	fg = rep(0, fm)
	fg[widx] = g
	if(optvars$index[2]>0) fg[midx] = -sum(2*optvars$benchmark*xport)/N
	return( fg )
}
##################################################################
# MAD
##################################################################
func.minmad = function(w, optvars, uservars)
{
	widx = optvars$widx
	f = mean( func.abs.smooth( optvars$Data %*% w[widx]- optvars$benchmark) )
	return( f )
}

grad.minmad = function(w, optvars, uservars)
{
	Data = optvars$Data
	widx = optvars$widx
	wm = optvars$wm
	fm = optvars$fm
	N = optvars$N
	g = vector(mode = "numeric", length = wm)
	xport = Data %*% w[widx] - optvars$benchmark
	for(i in 1:wm) g[i] = sum( Data[,i] * (parma.Small + xport)/(N * sqrt((parma.Small + xport)^2)))
	fg = rep(0, fm)
	fg[widx] = g
	return( fg )
}



func.optmad = function(w, optvars, uservars)
{	
	widx = optvars$widx
	midx = optvars$midx
	f =  mean( func.abs.smooth( optvars$Data %*% as.numeric(w[widx]) - w[midx] * optvars$benchmark) )
	return( f )
}


grad.optmad = function(w, optvars, uservars)
{
	Data = optvars$Data
	wm = optvars$wm
	fm = optvars$fm
	N = optvars$N
	widx = optvars$widx
	midx = optvars$midx
	g = rep(0, wm)
	xport = Data %*% w[widx] - w[midx]* optvars$benchmark
	for(i in 1:wm) g[i] = sum( Data[,i] * (parma.Small + xport)/(N * sqrt((parma.Small + xport)^2)))
	fg = rep(0, fm)
	fg[widx] = g
	if(optvars$index[2]>0) fg[midx] = -sum(optvars$benchmark * (parma.Small + xport)/(N * sqrt((parma.Small + xport)^2)))
	return( fg )
}
##################################################################
# CVaR
##################################################################

func.mincvar = function(w, optvars, uservars)
{
	Data = optvars$Data
	widx = optvars$widx
	vidx = optvars$vidx
	d = as.numeric( w[vidx] - (Data %*% w[widx]  - optvars$benchmark ))
	f = -w[vidx] + mean( func.max.smooth(d) )/optvars$alpha
	return( f )
}


grad.mincvar = function(w, optvars, uservars)
{
	Data = optvars$Data
	widx = optvars$widx
	vidx = optvars$vidx
	wm = optvars$wm
	fm = optvars$fm
	N = optvars$N
	alpha = optvars$alpha
	g = vector(mode = "numeric", length = wm)
	port = Data %*% w[widx]
	for(i in 1:wm){
		g[i] =  ( sum(Data[,i]/(2*N)) - sum( (Data[,i] * ( -optvars$benchmark + port - w[vidx])   )/( ( 2*N ) * 
										(parma.Small + ( -optvars$benchmark + port - w[vidx] )^2 )^0.5 ) ) )/alpha
	}
	fg = rep(0, fm)
	fg[widx] = g
	# derivative wrt to VaR
	fg[vidx] = ( sum( ( -optvars$benchmark + port - w[vidx] )/( 2*N * sqrt( parma.Small + ( -optvars$benchmark + port - w[vidx] )^2 ) ) ) - 0.5 )/alpha + 1 
	# we reverse the sign since we are minimizing
	return(-1 * fg)
}


func.optcvar = function(w, optvars, uservars){
	Data = optvars$Data
	widx = optvars$widx
	vidx = optvars$vidx
	midx = optvars$midx
	alpha = optvars$alpha	
	d = w[vidx] - (Data %*% w[widx] - w[midx]*optvars$benchmark)
	f = -w[vidx] + mean( func.max.smooth(d) )/alpha
	return( f )
}

grad.optcvar = function(w, optvars, uservars)
{
	Data = optvars$Data
	wm = optvars$wm
	fm = optvars$fm
	N = optvars$N
	widx = optvars$widx
	midx = optvars$midx
	vidx = optvars$vidx
	alpha = optvars$alpha
	g = rep(0, wm)
	port = Data %*% w[widx]
	for(i in 1:wm){
		g[i] =  ( sum(Data[,i]/(2*N)) - sum( (Data[,i] * ( port - w[vidx] - w[midx]*optvars$benchmark) )/( ( 2*N ) 
										* (parma.Small + ( port - w[vidx] - w[midx]*optvars$benchmark)^2 )^0.5 ) ) )/alpha
	}
	fg = rep(0, fm)
	fg[widx] = g
	# the derivative wrt to VaR
	fg[vidx] = ( sum( ( port - w[vidx] - w[midx]*optvars$benchmark)/( 2*N * sqrt( parma.Small 
													+ ( port - w[vidx] - w[midx]*optvars$benchmark)^2 ) ) ) - 0.5 )/alpha + 1
	if(optvars$index[2]>0) fg[midx] = sum(-optvars$benchmark/(2*N) + optvars$benchmark * 
						(port - w[vidx] - w[midx]*optvars$benchmark)/( ( 2*N ) * (parma.Small + ( port - w[vidx] - w[midx]*optvars$benchmark)^2 )^0.5))/(alpha)
	# we reverse the sign since we are minimizing
	return(-1 * fg)
}
##################################################################
# LPM
##################################################################
func.minlpm = function(w, optvars, uservars){
	moment = optvars$moment
	threshold = optvars$threshold
	widx = optvars$widx
	port = optvars$Data %*% as.numeric(w[widx])
	f =  (mean(func.max.smooth(threshold - port)^moment))^(1/moment)
	return( f )
}

grad.minlpm = function(w, optvars, uservars){
	Data = optvars$Data
	widx = optvars$widx
	wm = optvars$wm
	fm = optvars$fm
	N = optvars$N
	moment = optvars$moment
	threshold = optvars$threshold
	g = vector(mode = "numeric", length = wm)
	port = Data %*% as.numeric(w[widx])
	# avoid excessive calculations by precalculating some values:
	pd = port - threshold
	pd2 = pd*pd
	for(i in 1:wm){
		g[i] =  - ( ( sum( (abs( threshold + (parma.Small + pd2)^(1/2) - port))^(moment) )/(N*2^moment) )^(1/moment - 1) * 
					sum( moment * (Data[,i] - (Data[,i] * pd)/(parma.Small + pd2)^(1/2)) * 
									(abs(threshold + (parma.Small + pd2)^(1/2) - port))^(moment-1) ) )/(N*2^moment*moment)	
	}
	fg = rep(0, fm)
	fg[widx] = g
	return( fg )
}



func.optlpm = function(w, optvars, uservars){
	widx = optvars$widx
	midx = optvars$midx	
	moment = optvars$moment
	threshold = optvars$threshold
	port = optvars$Data %*% as.numeric(w[widx])
	f =  (mean(func.max.smooth( (w[midx]*threshold) - port)^moment))^(1/moment)
	return( f )
}

grad.optlpm = function(w, optvars, uservars){
	Data = optvars$Data
	wm = optvars$wm
	fm = optvars$fm
	N = optvars$N
	widx = optvars$widx
	midx = optvars$midx
	g = rep(0, wm)
	moment = optvars$moment
	threshold = optvars$threshold
	port = Data %*% as.numeric(w[widx])	
	# avoid excessive calculations by precalculating some values:
	tm = threshold*w[midx]
	pd = port - tm
	dp = tm - port
	pd2 = pd*pd	
	for(i in 1:wm){
		g[i] =  - (( sum( (abs(dp + sqrt(parma.Small + pd2)))^moment )/(N*2^moment) )^(1/moment - 1) *
					sum( moment * (Data[,i] - (Data[,i]*pd)/sqrt(parma.Small + pd2))*
									(abs(dp +sqrt(parma.Small + pd2)))^(moment-1)))/(moment*N*2^moment)
	}
	fg = rep(0, fm)
	fg[widx] = g
	# This is required because of the scaling of the threshold by the fractional parameter
	fg[midx] = (( sum( (dp + sqrt(parma.Small + (port-tm)^2))^moment )/(N*2^moment) )^(1/moment - 1) *
				sum(moment * (threshold - (threshold*pd)/sqrt(parma.Small + (port-tm)^2)) * 
								(dp + sqrt(parma.Small + pd2))^(moment-1)))/(N*moment*2^moment)
	return( fg )
}
##################################################################
# MiniMax
##################################################################
func.minminmax = function(w, optvars, uservars){
	return( w[optvars$vidx] )
}

grad.minminmax = function(w, optvars, uservars)
{
	fm = optvars$fm
	fg = rep(0, fm)
	fg[optvars$vidx] = 1
	return( fg )
}

func.optminmax = function(w, optvars, uservars){	
	return( w[optvars$vidx] )
}

grad.optminmax = function(w, optvars, uservars)
{
	fm = optvars$fm
	fg = rep(0, fm)
	fg[optvars$vidx] = 1
	return( fg )
}

ineqfun.minminmax = function(w, optvars, uservars)
{
	as.numeric( -w[optvars$vidx] - (as.numeric(optvars$Data %*% w[optvars$widx]) - optvars$benchmark) )
}

ineqjac.minminmax = function(w, optvars, uservars)
{
	g = cbind(matrix(-1, nrow = optvars$N), -optvars$Data)
	return(g)
}


ineqfun.optminmax = function(w, optvars, uservars)
{
	as.numeric( -w[optvars$vidx] - optvars$Data %*% w[optvars$widx] + w[optvars$midx]*optvars$benchmark )
}

ineqjac.optminmax = function(w, optvars, uservars)
{
	g = cbind(matrix(-1, nrow = optvars$N), -optvars$Data, matrix(-optvars$benchmark, nrow = optvars$N))
	return(g)
}
##################################################################
# LPM-UPM
##################################################################
func.optlpmupm = function(w, optvars, uservars){
	Data = optvars$Data
	widx = optvars$widx
	midx = optvars$midx
	lmoment = optvars$lmoment
	lthreshold = optvars$lthreshold
	port = Data %*% as.numeric(w[widx])
	f =  (mean(func.max.smooth( (w[midx]*lthreshold) - port)^lmoment))^(1/lmoment)
	return( f )
}
grad.optlpmupm = function(w, optvars, uservars){
	Data = optvars$Data
	wm = optvars$wm
	fm = optvars$fm
	N = optvars$N
	widx = optvars$widx
	midx = optvars$midx
	lmoment = optvars$lmoment
	lthreshold = optvars$lthreshold
	N = optvars$N
	g = rep(0, wm)
	port = Data %*% as.numeric(w[widx])	
	# avoid excessive calculations by precalculating some values:
	tm = lthreshold*w[midx]
	pd = port - tm
	dp = tm - port
	pd2 = pd*pd	
	for(i in 1:wm){
		g[i] =  - (( sum( (dp + sqrt(parma.Small + pd2))^lmoment )/(N*2^lmoment) )^(1/lmoment - 1) *
					sum( lmoment * (Data[,i] - (Data[,i]*pd)/sqrt(parma.Small + pd2))*
									(dp +sqrt(parma.Small + pd2))^(lmoment-1)))/(lmoment*N*2^lmoment)
	}
	fg = rep(0, fm)
	fg[widx] = g
	fg[midx] = (( sum( (dp + sqrt(parma.Small + (port-tm)^2))^lmoment )/(N*2^lmoment) )^(1/lmoment - 1) *
				sum(lmoment * (lthreshold - (lthreshold*pd)/sqrt(parma.Small + (port-tm)^2)) * 
								(dp + sqrt(parma.Small + pd2))^(lmoment-1)))/(N*lmoment*2^lmoment)
	return( fg )
}


func.ineq.lpmupm = function(w, optvars, uservars){
	Data = optvars$Data
	widx = optvars$widx
	midx = optvars$midx
	umoment = optvars$umoment
	uthreshold = optvars$uthreshold
	N = optvars$N
	LB = optvars$LB
	UB = optvars$UB
	port = Data %*% as.numeric(w[widx])
	A1 =  1 - (mean(func.max.smooth(port - uthreshold*w[midx])^umoment))^(1/umoment)
	A2 = ( c( w[midx]*LB - w[widx], w[widx] - w[midx]*UB) )
	return( c(A1, A2) )
}


grad.ineq.lpmupm = function(w, optvars, uservars){
	Data = optvars$Data
	widx = optvars$widx
	midx = optvars$midx
	wm = optvars$wm
	fm = optvars$fm
	umoment = optvars$umoment
	uthreshold = optvars$uthreshold
	N = optvars$N
	LB = optvars$LB
	UB = optvars$UB
	g = rep(0, wm)
	port = Data %*% as.numeric(w[widx])
	# avoid excessive calculations by precalculating some values:
	pd =  port - uthreshold*w[midx]
	pd2 = pd*pd
	for(i in 1:wm){
		ans1[i] =  -(( sum( (pd + sqrt(parma.Small + pd2))^umoment )/(N*2^umoment) )^(1/umoment - 1) *
					sum( umoment * (Data[,i] + (Data[,i]*pd)/sqrt(parma.Small + pd2))*(pd +sqrt(parma.Small + pd2))^(umoment-1)))/(umoment*N*2^umoment)
	}
	fg = rep(0, fm)
	fg[widx] = g
	
	fg[midx] = (  ( (sum( (pd + sqrt(parma.Small + pd2))^umoment)/(N*2^umoment))^(1/umoment-1) ) * 
				(  sum( umoment*(uthreshold + (uthreshold*pd)/sqrt(parma.Small+pd2))*(pd + sqrt(parma.Small + pd2))^(umoment-1))))/(umoment*N*2^umoment)
	
	xl = xu = matrix(0, nrow = fm, ncol = fm)
	xu[widx, widx] = -diag(wm)
	xl[widx, widx] =  diag(wm)
	xu[widx, midx] = -UB	
	xl[widx, midx] =  LB
	xu = xu[widx, ]
	xl = xl[widx, ]
	
	ans = rbind(matrix(fg, nrow = 1), xl, xu)
	return(ans)
}

##################################################################
# CONSTRAINTS
##################################################################
#-----------------------------------------------------
# Budget
#-----------------------------------------------------
eqfun.budget.min = function(w, optvars, uservars){
	widx = optvars$widx	
	return( sum(w[widx]) - optvars$budget )
}

eqjac.budget.min = function(w, optvars, uservars){
	widx = optvars$widx
	fm = optvars$fm
	j = matrix(0, 1, fm)
	j[1, widx] = 1
	return(j)
}

eqfun.budget.opt = function(w, optvars, uservars){
	midx = optvars$midx
	widx = optvars$widx
	return( sum(w[widx]) - w[midx]*optvars$budget )
}

eqjac.budget.opt = function(w, optvars, uservars){
	widx = optvars$widx
	midx = optvars$midx
	fm = optvars$fm
	j = matrix(0, 1, fm)
	j[1, widx] = 1
	j[1, midx] = -optvars$budget
	return(j)
}
#-----------------------------------------------------
# Leverage
#-----------------------------------------------------
eqfun.leverage.min = function(w, optvars, uservars){
	widx = optvars$widx	
	return( sum( sqrt(w[widx]^2 + parma.Small) ) - optvars$leverage )
}

eqjac.leverage.min = function(w, optvars, uservars){
	widx = optvars$widx
	fm = optvars$fm
	j = matrix(0, 1, fm)
	j[1, widx] = w[widx]/sqrt(w[widx]^2 + parma.Small)
	return(j)
}

eqfun.leverage.opt = function(w, optvars, uservars){
	widx = optvars$widx
	midx = optvars$midx
	return( sum( sqrt(w[widx]^2 + parma.Small) ) - w[midx]*optvars$leverage )
}

eqjac.leverage.opt = function(w, optvars, uservars){
	widx = optvars$widx
	midx = optvars$midx
	fm = optvars$fm
	j = matrix(0, 1, fm)
	j[1, widx] = w[widx]/sqrt(w[widx]^2 + parma.Small)
	j[1, midx] = -optvars$leverage
	return(j)
}

#-----------------------------------------------------
# Target
#-----------------------------------------------------
# Equality
# if benchmark relative optimization:-> 
# mu = excess return and target = excess target
eqfun.target.min = function(w, optvars, uservars){
	widx = optvars$widx	
	return( sum(w[widx] * optvars$mu) - optvars$mutarget)
}

eqjac.target.min = function(w, optvars, uservars){
	widx = optvars$widx
	fm = optvars$fm
	j = matrix(0, 1, fm)
	j[1, widx] = optvars$mu
	return(j)	
}
eqfun.target.opt = function(w, optvars, uservars){
	widx = optvars$widx
	midx = optvars$midx
	return( sum(w[widx]*optvars$mu) - 1 )
}

eqjac.target.opt = function(w, optvars, uservars){
	widx = optvars$widx
	midx = optvars$midx
	fm = optvars$fm
	j = matrix(0, 1, fm)
	j[1, widx] = optvars$mu
	#j[1, midx] = -optvars$mbenchmark
	return(j)
}

# Inequality g(w)<=0 (only for minrisk since optimal is not valid)
ineqfun.target.min = function(w, optvars, uservars){
	widx = optvars$widx	
	return( optvars$mutarget - sum(w[widx] * optvars$mu) )
}

ineqjac.target.min = function(w, optvars, uservars){
	widx = optvars$widx
	fm = optvars$fm
	j = matrix(0, 1, fm)
	j[1, widx] = -optvars$mu
	return(j)
}
#-----------------------------------------------------
# Bounds (fractional problem)
#-----------------------------------------------------
ineqfun.bounds.opt = function(w, optvars, uservars){
	midx = optvars$midx
	widx = optvars$widx
	LB = optvars$LB
	UB = optvars$UB
	return( c( w[midx]*LB - w[widx], w[widx] - w[midx]*UB ) )
}

ineqjac.bounds.opt = function(w, optvars, uservars){
	midx = optvars$midx
	widx = optvars$widx
	fm = optvars$fm
	LB = optvars$LB
	UB = optvars$UB
	m = length(widx)
	ans = matrix(0, nrow = 2*m, ncol = fm)
	ans[1:m, widx] = -diag(m)
	ans[(m+1):(2*m), widx] = diag(m)
	ans[1:m, midx] = LB
	ans[(m+1):(2*m), midx] = -UB
	return( ans )
}

#-----------------------------------------------------
# LPM Inequality (LPM-UPM Problem)
#-----------------------------------------------------
ineqfun.lpm.min = function(w, optvars, uservars){
	Data = optvars$Data
	widx = optvars$widx
	moment = uservars$lmoment
	threshold = uservars$lthreshold
	port = Data %*% as.numeric(w[widx])
	f =  (mean(func.max.smooth( threshold - port)^moment))^(1/moment)
	return(f - uservars$lpmtarget)
}


ineqjac.lpm.min = function(w, optvars, uservars){
	Data = optvars$Data
	widx = optvars$widx
	wm = optvars$wm
	fm = optvars$fm
	N = optvars$N
	moment = uservars$lmoment
	threshold = uservars$lthreshold
	g = vector(mode = "numeric", length = wm)
	port = Data %*% as.numeric(w[widx])
	# avoid excessive calculations by precalculating some values:
	pd = port - threshold
	pd2 = pd*pd
	for(i in 1:wm){
		g[i] =  - ( ( sum( ( threshold + (parma.Small + pd2)^(1/2) - port)^(moment) )/(N*2^moment) )^(1/moment - 1) * 
					sum( moment * (Data[,i] - (Data[,i] * pd)/(parma.Small + pd2)^(1/2)) * 
									(threshold + (parma.Small + pd2)^(1/2) - port)^(moment-1) ) )/(N*2^moment*moment)	
	}
	fg = rep(0, fm)
	fg[widx] = g
	return( matrix(fg, nrow = 1) )
}


ineqfun.lpm.opt = function(w, optvars, uservars){
	Data = optvars$Data
	widx = optvars$widx
	midx = optvars$midx
	moment = uservars$lmoment
	threshold = uservars$lthreshold
	port = Data %*% as.numeric(w[widx])
	f =  (mean(func.max.smooth( w[midx]*threshold - port)^moment))^(1/moment)
	return(f - w[midx]*uservars$lpmtarget)
}



ineqjac.lpm.opt = function(w, optvars, uservars){
	Data = optvars$Data
	wm = optvars$wm
	fm = optvars$fm
	N = optvars$N
	widx = optvars$widx
	midx = optvars$midx
	g = rep(0, wm)
	moment = uservars$lmoment
	threshold = uservars$lthreshold
	port = Data %*% as.numeric(w[widx])	
	# avoid excessive calculations by precalculating some values:
	tm = threshold*w[midx]
	pd = port - tm
	dp = tm - port
	pd2 = pd*pd	
	for(i in 1:wm){
		g[i] =  - (( sum( (dp + sqrt(parma.Small + pd2))^moment )/(N*2^moment) )^(1/moment - 1) *
					sum( moment * (Data[,i] - (Data[,i]*pd)/sqrt(parma.Small + pd2))*
									(dp +sqrt(parma.Small + pd2))^(moment-1)))/(moment*N*2^moment)
	}
	fg = rep(0, fm)
	fg[widx] = g
	# This is required because of the scaling of the threshold by the fractional parameter
	fg[midx] = -uservars$lpmtarget + (( sum( (dp + sqrt(parma.Small + (port-tm)^2))^moment )/(N*2^moment) )^(1/moment - 1) *
				sum(moment * (threshold - (threshold*pd)/sqrt(parma.Small + (port-tm)^2)) * 
								(dp + sqrt(parma.Small + pd2))^(moment-1)))/(N*moment*2^moment)
	
	return( matrix(fg, nrow = 1) )
}
################################################################
nlpport = function(optvars, uservars, control = list(), ...)
{
	ctrl = .nloptr.ctrl(control)
	xrisk = c("mad", "minimax", "cvar", "cdar", "ev", "lpm", "lpmupm")[optvars$index[4]]
	sol = list()

	if(optvars$index[5] == 1){
		if(optvars$index[3]==2 && is.null(optvars$ineqfun) && xrisk!="minimax"){
			xparma.ineq.minfun = NULL
			xparma.ineq.mingrad = NULL
		} else{
			xparma.ineq.minfun = parma.ineq.minfun
			xparma.ineq.mingrad = parma.ineq.mingrad
		}
		eval_f = switch(xrisk, 
				"cvar" = func.mincvar,
				"mad"  = func.minmad, 
				"lpm"  = func.minlpm,
				"minimax" = func.minminmax,
				"ev" = func.minvar)
		
		eval_grad_f = switch(xrisk, 
				"cvar" = grad.mincvar,
				"mad"  = grad.minmad, 
				"lpm"  = grad.minlpm,
				"minimax" = grad.minminmax,
				"ev" = grad.minvar)
		
		ans = try(nloptr(x0 = optvars$x0, 
				eval_f = eval_f, 
				eval_grad_f = eval_grad_f, 
				lb = optvars$LB, ub = optvars$UB, 
				eval_g_ineq = xparma.ineq.minfun, 
				eval_jac_g_ineq = xparma.ineq.mingrad, 
				eval_g_eq = parma.eq.minfun, 
				eval_jac_g_eq = parma.eq.mingrad, 
				opts = ctrl, 
				optvars = optvars, uservars = uservars), silent = TRUE)
		
		if(inherits(ans, 'try-error')){
			print(ans)
			sol$status = 1
			sol$weights = rep(NA, optvars$wm)
			sol$risk = NA
			sol$reward = NA
			if(optvars$index[4]==3) sol$VaR = NA
			if(optvars$index[4]==4) sol$DaR = NA
		} else{
			sol$status = ans$status
			sol$weights = ans$solution[optvars$widx]
			sol$risk = as.numeric(fun.risk(sol$weights, optvars$Data, 
							options = list(alpha = optvars$alpha, 
									moment = optvars$moment, 
									threshold = optvars$threshold), 
							risk = xrisk, benchmark = optvars$benchmark))
			sol$reward = sum(sol$weights * optvars$mu)
			if(optvars$index[4]==3) sol$VaR = ans$solution[optvars$vidx]
			if(optvars$index[4]==4) sol$DaR = ans$solution[optvars$vidx]
		}
	} else{
		eval_f = switch(xrisk, 
				"cvar" = func.optcvar,
				"mad"  = func.optmad, 
				"lpm"  = func.optlpm,
				"minimax" = func.optminmax,
				"ev" = func.optvar)
		
		eval_grad_f = switch(xrisk, 
				"cvar" = grad.optcvar,
				"mad"  = grad.optmad, 
				"lpm"  = grad.optlpm,
				"minimax" = grad.optminmax,
				"ev" = grad.optvar)
		
		ans = try(nloptr(x0 = optvars$x0, 
				eval_f = eval_f, 
				eval_grad_f = eval_grad_f, 
				lb = optvars$fLB, ub = optvars$fUB, 
				eval_g_ineq = parma.ineq.optfun, 
				eval_jac_g_ineq = parma.ineq.optgrad, 
				eval_g_eq = parma.eq.optfun, 
				eval_jac_g_eq = parma.eq.optgrad, 
				opts = ctrl, 
				optvars = optvars, uservars = uservars), silent = TRUE)

		if(inherits(ans, 'try-error')){
			print(ans)
			sol$status = 1
			sol$weights = rep(NA, optvars$wm)
			sol$risk = NA
			sol$reward = NA
			sol$multiplier = NA
			if(optvars$index[4]==3) sol$VaR = NA
			if(optvars$index[4]==4) sol$DaR = NA
		} else{
			sol$status = ans$status
			sol$weights = ans$solution[optvars$widx]/ans$solution[optvars$midx]
			sol$reward = sum(sol$weights * optvars$mu)
			sol$multiplier = ans$solution[optvars$midx]
			sol$risk = as.numeric(fun.risk(sol$weights, optvars$Data, options = list(alpha = optvars$alpha, 
							moment = optvars$moment, threshold = optvars$threshold), 
					risk = xrisk, benchmark = optvars$benchmark))
			if(optvars$index[4]==3) sol$VaR = ans$solution[optvars$vidx]/ans$solution[optvars$midx]
			if(optvars$index[4]==4) sol$DaR = ans$solution[optvars$vidx]/ans$solution[optvars$midx]
		}
	}
	return( sol )
}

.nloptr.ctrl = function(control = list()){
	ctrl = list()
	ctrl$algorithm = "NLOPT_LD_SLSQP"
	if( is.null(control$minf_max) ) ctrl$minf_max = 1e-24 else ctrl$minf_max = control$minf_max
	if( is.null(control$ftol_rel) ) ctrl$ftol_rel = 1e-16 else ctrl$ftol_rel = control$ftol_rel
	if( is.null(control$ftol_abs) ) ctrl$ftol_abs = 1e-16 else  ctrl$ftol_abs = control$ftol_abs
	if( is.null(control$xtol_rel) ) ctrl$xtol_rel = 1e-24 else  ctrl$xtol_rel = control$xtol_rel
	if( is.null(control$maxeval) ) ctrl$maxeval = 10000 else  ctrl$maxeval = as.integer( control$maxeval )
	if( is.null(control$maxtime) ) ctrl$maxtime = 10000 else  ctrl$maxtime = control$maxtime
	if( is.null(control$print_level) ) ctrl$print_level = 0 else  ctrl$print_level = as.integer( control$print_level )
	return(ctrl)
}

###############################################################################
# Some Custom Constraints:

#-----------------------------------------------------
# Turnover
#-----------------------------------------------------
# Minimum Risk
ineqfun.turnover.min = function(w, optvars, uservars){
	widx = optvars$widx
	wold = uservars$wold
	return( sum( func.abs.smooth(w[widx] - wold) ) - uservars$turnover )
}

ineqjac.turnover.min = function(w, optvars, uservars){
	widx = optvars$widx
	fm = optvars$fm
	j = matrix(0, 1, fm)
	j[1, widx] = grad.abs.smooth(w[widx]- uservars$wold)
	return(j)
}

ineqfun.bsturnover.min = function(w, optvars, uservars){
	widx = optvars$widx
	wold = uservars$wold
	bt = sum(func.max.smooth(w[widx] - wold)) - uservars$buyturnover
	st = sum(func.max.smooth(wold - w[widx])) - uservars$sellturnover
	return(c(bt, st))
}

ineqjac.bsturnover.min = function(w, optvars, uservars){
	widx = optvars$widx
	fm = optvars$fm
	j = matrix(0, 2, fm)
	j[1, widx] = 0.25*(-2*uservars$wold+2*w[widx])/sqrt( (-uservars$wold + w[widx])^2 + parma.Small ) + 0.5
	j[2, widx] = 0.25*(-2*uservars$wold+2*w[widx])/sqrt( (uservars$wold - w[widx])^2 + parma.Small ) - 0.5
	return(j)
}

#-------------------------------------------------------------------------------
# Fractional
ineqfun.turnover.opt = function(w, optvars, uservars){
	widx = optvars$widx
	midx = optvars$midx
	wold = uservars$wold
	print(sum(w[widx]/w[midx]))
	return( sum( func.abs.smooth(w[widx]/w[midx] - wold) ) - uservars$turnover )
}

ineqjac.turnover.opt = function(w, optvars, uservars){
	widx = optvars$widx
	midx = optvars$midx
	fm = optvars$fm
	j = matrix(0, 1, fm)
	tmp = w[widx]/w[midx] - uservars$wold + 1e-20
	j[1, widx] = tmp/(w[midx]*sqrt(tmp^2))
	j[1, midx] = sum(- ( (tmp*w[widx])/( (w[midx]^2)*sqrt(tmp^2) ) ) )
	return(j)
}


ineqfun.bsturnover.opt = function(w, optvars, uservars){
	widx = optvars$widx
	wold = uservars$wold
	midx = optvars$midx
	bt = sum(func.max.smooth(w[widx]/w[midx] - wold)) - uservars$buyturnover
	st = sum(func.max.smooth(wold - w[widx]/w[midx])) - uservars$sellturnover
	return(c(bt, st))
}

ineqjac.bsturnover.opt = function(w, optvars, uservars){
	widx = optvars$widx
	midx = optvars$midx
	fm = optvars$fm
	j = matrix(0, 2, fm)
	tmp = w[widx]/w[midx] - uservars$wold
	j[1, widx] = 0.5*(tmp/(w[midx]*sqrt(tmp^2+1e-20))) + 0.5/w[midx]
	j[1, midx] = sum( - 0.5*tmp*w[widx]/(w[midx]^2*sqrt(tmp^2+1e-20)) - 0.5*w[widx]/w[midx]^2 )
	tmp = -w[widx]/w[midx] + uservars$wold
	j[2, widx] = -0.5*(tmp/(w[midx]*sqrt(tmp^2+1e-20))) - 0.5/w[midx]
	j[2, midx] = sum( 0.5*tmp*w[widx]/(w[midx]^2*sqrt(tmp^2+1e-20)) + 0.5*w[widx]/w[midx]^2 )
	return(j)
}
ineqfun.variance.opt = function(w, optvars, uservars)
{
	# rescaled weights
	x = w[optvars$widx]/w[optvars$midx]
	return( as.numeric( x %*% uservars$Cov %*% x) - uservars$varbound )
}
ineqjac.variance.opt = function(w, optvars, uservars){
	g = matrix(0, ncol = optvars$fm, nrow = 1)
	widx = optvars$widx
	midx = optvars$midx
	x = w[widx]
	tmp1 = w[widx]/w[midx]
	tmp2 = w[widx]/w[midx]^2
	g[1, widx] = as.numeric( ( tmp1 %*% uservars$Cov)/w[midx] + (x %*% uservars$Cov)/w[midx]^2 )
	g[1, midx] = as.numeric( ( (tmp2/w[midx]) %*% -uservars$Cov)%*%w[widx]- ( (tmp1/w[midx]^2) %*% uservars$Cov )%*%w[widx] )
	return(g)
}


ineqfun.variance.min = function(w, optvars, uservars)
{
	x = w[optvars$widx]
	return( as.numeric( x %*% uservars$Cov %*% x) - uservars$varbound )
}
ineqjac.variance.min = function(w, optvars, uservars){
	g = matrix(0, ncol = optvars$fm, nrow = 1)
	widx = optvars$widx
	midx = optvars$midx
	x = w[widx]
	g[1, widx] = as.numeric( 2*x %*% uservars$Cov )
	return(g)
}
