#################################################################################
##
##   R package spd by Alexios Ghalanos Copyright (C) 2008-2013
##   This file is part of the R package spd.
##
##   The R package spd is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package spd is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################
# approximate PDF from CDF given a set of bin points found by calling
# .histcount
eps=.Machine$double.eps
.interpPDF = function(xm,ppx,ppz)
{
	# xm = middle values of dspd(x,...)
	# ppx = kernel interpolated density of original data
	# ppz = empirical cdf probability [0,1] (cumsum of ppy/sum of ppy)
	m=length(xm)
	bins<-vector(mode="numeric",length=m)
	nbin = length(ppx)
	nx = vector(mode="numeric",length=nbin)
	z=as.data.frame(1:nbin)
	kk=apply(z,1,FUN=function(x) which( xm >= ppx[x]))
	for(i in 1:nbin)
	{
		bins[kk[[i]]]<-i
	}
	nx=apply(z,1,FUN=function(x) length(kk[[x]]))
	kk = which( xm > ppx[nbin] )
	bins[kk] = 0
	nx[nbin+1] = length(kk)
	counts = -diff(nx)
	bin=bins
	bin[xm==ppx[1]] = 1
	bin[xm==ppx[length(ppx)]] = length(counts)-1
	fx = vector(mode="numeric",length=length(xm))
	tx = bin>0
	bin = bin[tx]
	# find the PDF probability
	tmp = (ppz[bin+1] - ppz[bin])/ (ppx[bin+1] - ppx[bin])
	tna = which(!is.na(tmp))
	tx = tx[tna]
	fx[tna] = tmp[!is.na(tmp)]
	return(fx)
}

.findthresh<-function(data, exceed)
{
	data <- rev(sort(data))
	uniq <- unique(data)
	idx <- match(data[exceed], uniq)
	idx <- pmin(idx + 1, length(uniq))
	return(uniq[idx])
}

.description<-function() 
{
	ans = paste(as.character(date()), "by user:", Sys.getenv("USERNAME"))
	ans
}
