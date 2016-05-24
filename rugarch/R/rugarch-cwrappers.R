#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2015.
##   This file is part of the R package rugarch.
##
##   The R package rugarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rugarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################
# fractional difference series C wrapper
.arfimaxfilter = function(model, pars, idx, mexdata, h, data, N)
{
	#if(model[1] == 0) pars[1,1] = 0
	m = as.integer(N[1])
	T = as.integer(N[2])
	if(length(h) <= 1) {
		h = double(length = T)
	} else{
		h = as.double(h)
	}
	data = as.double(data)
	# flatten exogenous matrix
	if(model[6]>0){
		xmxreg = matrix( pars[idx[6,1]:idx[6,2]], ncol = model[6] )
		if(model[20]==0){
			imx =  xmxreg %*%t( matrix( mexdata, ncol = model[6] ) )
		} else{
			if(model[20] == model[6]){
				imx = xmxreg %*%t( matrix( mexdata * h , ncol = model[6] ) )				
			} else{
				imx = xmxreg[,1:(model[6]-model[20]),drop=FALSE] %*%t( matrix( mexdata[,1:(model[6]-model[20]),drop=FALSE], ncol = model[6]-model[20] ) )
				imx = imx + xmxreg[,(model[6]-model[20]+1):model[6],drop=FALSE] %*%t( matrix( mexdata[,(model[6]-model[20]+1):model[6],drop=FALSE]*h, ncol = model[20] ) )					
			}
		}
		imx = as.double(imx)
		mexdata = as.double(as.vector(mexdata))
	} else{
		mexdata = as.double(0)
		imx = as.double(0)
	}
	res = double(length = T)
	# this routine is used for the mean residuals to initiate the recursion
	# so we ignore arfima before
	zrf = double(length = T)
	constm = double(length = T)
	condm = double(length = T)
	ans = list()
	if(model[2]>0 | model[3]>0){
		ans = try(.C("arfimaxfilterC", model = as.integer(model), pars = as.double(pars), 
						idx = as.integer(idx[,1]-1), x = data, res = res, mexdata = mexdata, 
						zrf = zrf, constm = constm, condm = condm, h = h, m = m, T = T,  
						PACKAGE = "rugarch"), silent = TRUE)
		if(inherits(ans, "try-error") | any(is.nan(ans$res)) | any(is.na(ans$res)) | any(!is.finite(ans$res)) ){
			res = data - pars[idx[1,1]]
			ans$res = res
			if(model[4]>0)
			{
				ans$zrf = .fracdiff(c(1,rep(0,length(data)-1)), darfima = pars[idx[4,1]])
				ans$res = .fracdiff(ans$res, darfima = pars[idx[4,1]])
			}
			if(any(is.na(res))) res[which(is.na(res))]=0
			return(ans)
		} else{
			if(model[4]>0)
			{
				ans$zrf = .fracdiff(c(1, rep(0,length(data)-1)), darfima = pars[idx[4,1]])
				ans$res = .fracdiff(ans$res, darfima = pars[idx[4,1]])
			}
			if(any(is.na(ans$res))) res[which(is.na(ans$res))]=0
			return(ans)
		}
	} else{
		ans = list()
		ans$res = data -  pars[idx[1,1]] - imx - pars[idx[5,1]]*(h^model[5])
		ans$zrf = zrf
		if(model[4]>0)
		{
			ans$zrf = .fracdiff(c(1,rep(0,length(data)-1)), darfima = pars[idx[4,1]])
			ans$res = .fracdiff(ans$res, darfima = pars[idx[4,1]])
		}
		if(any(is.na(ans$res))) res[which(is.na(ans$res))]=0
		return(ans)
	}
}


# fractional difference series C wrapper
.fracdiff = function(x, darfima)
{
	n = length(as.vector(x))
	p = c(-darfima, rep(0,n-1))
	
	res = .C("fracdiff",n = as.integer(n), d = as.double(darfima), p = as.double(p),
			x = as.double(x), ydiff = as.double(x), PACKAGE = "rugarch")
	return(res$ydiff)
}

.arfimafitC = function(model, pars, idx, mexdata, sigma, data, zrf, N, res)
{
	m = as.integer(N[1])
	T = as.integer(N[2])
	sigma = as.double(sigma)
	zrf = as.double(zrf)
	data = as.double(data)	
	# flatten exogenous matrix	
	if(model[6]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = as.double(0)	
	z = double(T)
	h = double(T)
	res = as.double(res)
	constm = double(T)
	condm = double(T)
	llh = double(1)
	LHT = double(T)
	
	ans = try(.C("arfimafitC", model = as.integer(model), pars = as.double(pars), 
					idx = as.integer(idx), sigma = sigma, x = data, res = res, 
					mexdata = mexdata, zrf = zrf, constm = constm, condm = condm, 
					m = m, T = T, h = h, z = z, llh = llh, LHT = LHT, 
					PACKAGE = "rugarch"), silent = TRUE)
	if(inherits(ans, "try-error")){
		return(0)
	} else{
		return(ans)
	}
}


.armaxsim = function(model, ipars, idx, constm, x, res, T, m)
{
	ans = try(.C("armaxsim", model = as.integer(model), pars = as.double(ipars[,1]), 
					idx = as.integer(idx[,1]-1), x = as.double(x), res = as.double(res), 
					constm = as.double(constm), m = as.integer(m), T = as.integer(T), 
					PACKAGE = "rugarch"), silent = TRUE)
	if(inherits(ans, "try-error")){
		return(0)
	} else{
		return(ans)
	}
}
.arfimaxsim = function(model, ipars, idx, constm, res, T)
{
	res = as.double(res)
	T = as.integer(T)
	constm = as.double(constm)
	flmin = as.double(.Machine$double.xmin)
	flmax = as.double(.Machine$double.xmax)
	epmin = as.double(.Machine$double.neg.eps)
	epmax = as.double(.Machine$double.eps)
	s = double(T+model[3])
	d = as.double( ipars[idx["arfima",1], 1] )
	d = max(1e-9, d)
	d = min(0.5-1e-9, d)
	ans = list()
	ans = try(.Fortran("fdsim", n = T, ip = as.integer( model[2] ), iq = as.integer( model[3] ), 
					ar = as.double( ipars[idx["ar",1]:idx["ar",2], 1] ), 
					ma = as.double( ipars[idx["ma",1]:idx["ma",2], 1] ), 
					d = as.double(d),
					rmu = constm, y = res, s = s, flmin = flmin, flmax = flmax,
					epmin = epmin, epmax = epmax,
			PACKAGE = "rugarch"), silent = TRUE)
	if(inherits(ans, "try-error")){
		return(0)
	} else{
		ans$series = ans$s
		ans$s = NULL
		return(ans)
	}
}