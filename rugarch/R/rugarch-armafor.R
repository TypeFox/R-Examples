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
# forecast from an arma/garch-in-mean process
armaf = function(ipars, model, idx, mu, mxfi, h, epsx, z, data, N, n.ahead)
{
	# mxf is the extended external regressor matrix to
	# include their n.ahead forecasts provided in the
	# forecasting routine stage:
	# mxf: rbind(mexdata, external.forecasts$mregfor)
	
	# 20-09-2011 bug report was already passing data+n.ahead in some routines
	# x = c(data, rep(0, n.ahead))
	x = data
	if(model[6]>0){
		if(model[20]==0){
			mu = mu + mxfi%*%t(matrix(ipars[idx["mxreg",1]:idx["mxreg",2],1], ncol = model[6]))
		} else{
			if(model[20] == model[6]){
				mu = mu + (mxfi*h)%*%t(matrix(ipars[idx["mxreg",1]:idx["mxreg",2],1], ncol = model[6]))
				
			} else{
				mu = mu + (mxfi[,1:(model[6]-model[20])])%*%t(matrix(ipars[idx["mxreg",1]:(idx["mxreg",1]+(model[6]-model[20]-1)),1], ncol = model[6] - model[20]))
				mu = mu + (mxfi[,(model[6]-model[20]+1):model[6]]*h)%*%t(matrix(tail(ipars[idx["mxreg",1]:idx["mxreg",2],1], model[20]), ncol = model[20]))
			}
		}
	}
	# comment: we allow the arch-in-mean to decay (towards its unconditional value) 
	# at the rate of the variance forecast (->1 - persistence) rather than just using 
	# E[h]*inmean...
	# also note that h here is sigma not the variance
	if(model[5]>0) mu = mu + ipars[idx["archm",1],1]*(h^model[5])
	# correction for presence of xreg ( 20 Apr.2011 )
	for(i in 1:n.ahead){
		x[N+i] = mu[N+i] + ifelse(model[2]>0, sum(ipars[idx["ar",1]:idx["ar",2],1]*(x[N+i-(1:model[2])] - mu[N+i-(1:model[2])])),0)
		if(model[3]>0){
			for(j in 1:model[3]){
				if(i-j>0){
					s = 0
				} else{
					s = ipars[idx["ma",1]+j-1,1]*epsx[N+i-j]
				}
				x[N+i] = x[N+i] + s
			}
		}
	}
	return(x)
}

# forecast with arfima
arfimaf = function(ipars, model, idx, mu, mxfi, h, epsx, z, data, N, n.ahead)
{
	if(model[6]>0){
		if(model[20]==0){
			mu = mu + mxfi%*%t(matrix(ipars[idx["mxreg",1]:idx["mxreg",2],1], ncol = model[6]))
		} else{
			if(model[20] == model[6]){
				mu = mu + (mxfi*h)%*%t(matrix(ipars[idx["mxreg",1]:idx["mxreg",2],1], ncol = model[6]))
				
			} else{
				mu = mu + (mxfi[,1:(model[6]-model[20])])%*%t(matrix(ipars[idx["mxreg",1]:(idx["mxreg",1]+(model[6]-model[20])),1], ncol = model[6] - model[20]))
				mu = mu + (mxfi[,(model[6]-model[20]+1):model[6]]*h)%*%t(matrix(tail(ipars[idx["mxreg",1]:idx["mxreg",2],1], model[20]), ncol = model[20]))
			}
		}
	}
	if(model[5]>0) mu = mu + ipars[idx["archm",1],1]*(h^model[5])
	# 20-09-2011 bug report was already passing data+n.ahead in some routines
	#x = c(data, rep(0, n.ahead))
	x = data
	# demean data (includes regressors and arch-in-mean)
	xdata = x - mu
	# N should now be as passed, not N+n.ahead=length(data)
	#N = length(data)
	# generate fractional series
	zrf = .fracdiff(c( 1, rep(0, N + n.ahead - 1) ), darfima = as.double( ipars[idx["arfima",1],1]) )
	if(model[2]>0){
		zar = rep(0, model[2])
		for(i in 1:model[2]){
			zar[i] = -sum(rev(zrf[2:(N - i + 1)])*xdata[1:(N-i)])
		}
		# reverse back
		zar = rev(zar)
		# lagged data adjusted for fractional integration (i.e. subtracted
		# fractional mean)
		xzr = xdata[1:N] - c(rep(0, N-model[2]), zar)
	} else{
		# no need in the absence of armap since we have no lagged x terms
		# i.e. no 0 * x-mu = 0
		xzr = xdata[1:N]
	}
	for(i in 1:n.ahead){
		# add fractional mean back
		xzy = -sum(rev(zrf[2:(N + i)])*xdata[1:(N+i-1)])
		x[N+i] = (mu[N+i] + xzy) + ifelse(model[2]>0, sum(ipars[idx["ar",1]:idx["ar",2],1]*(xzr[N+i-(1:model[2])])), 0)
		if(model[3]>0){
			for(j in 1:model[3]){
				if(i-j>0){
					s = 0
				} else{
					s = ipars[idx["ma",1]+j-1,1]*epsx[N+i-j]
				}
				x[N+i] = x[N+i] + s
			}
		}		
		# forward adjust for the mean
		xdata[N+i] = x[N+i] - mu[N+i]
		# forward adjust for the lagged data (irrelevant when armap=0)
		xzr[N+i] = x[N+i] - mu[N+i] - xzy
	}
	return(x)
}