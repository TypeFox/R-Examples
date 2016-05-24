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
# Developer Note: extendable class and method for different tail distributions.
# perhaps as we add more tailfit models it would make sense to create a class and
# method for creating a distribution design object which would hold the specific details
# which would then be passed on to the fit (i.e. distributionSpec)

spdfit<-function(data, upper = 0.9, lower = 0.1, tailfit="GPD", type = c("mle", "pwm"), kernelfit = c("normal","box","epanech","biweight","triweight"), information = c("observed", "expected"), title = NULL, description = NULL,...)
{
	UseMethod("spdfit")
}

.spdfit<-function(data, upper = 0.9, lower = 0.1, tailfit = "GPD", type = c("mle", "pwm"), kernelfit = c("normal","box","epanech","biweight","triweight"), information = c("observed", "expected"), title = NULL, description = NULL,...)
{
	ans<-switch(tailfit,
			GPD = .gpdtails(data, upper, lower, tailfit, type, kernelfit, information, title, description,...))
			#GEV = .gevtails(data, upper, lower, tailfit, type, kernelfit, information, title, description,...))
	return(ans)
}

setMethod(f="spdfit", definition=.spdfit)
#-------------------------------------------------------------------------------------------------------------------
.gpdtails<-function(data, upper = 0.9, lower = 0.1, tailfit="GPD", type = c("mle", "pwm"), kernelfit = c("normal","box","epanech","biweight","triweight"), information = c("observed", "expected"), title = NULL, description = NULL,...)
{
	# need to add validity checks for type, kernel and information
	if(!missing(kernelfit) & length(kernelfit)>1) kernelfit<-"normal"
	if(missing(kernelfit)) kernelfit<-"normal"
	if(is.null(title)) title=""
	if(is.null(description)) description=""
    x=data
    x = sort(as.matrix(x)[,1])
    call = match.call()
    type = match.arg(type)
    kernel=match.arg(kernelfit)
    information = match.arg(information)
    # Check Type and Convert:
    x = as.vector(x)
    N = length(x)
    lp = trunc(N*lower)
    up = trunc(N*(1-upper))
    nUpper = x[N - up]
    nLower = x[lp]
    #upper tail
    upperExceedances = x[x > nUpper]
    upperExcess = upperExceedances - nUpper
    #lower tail
    lowerExceedances = x[x < nLower]
    lowerExcess = nLower - lowerExceedances
    #estimate GPD Tails and kernel interior
    upperfit = gpdfit(x, nUpper, type = type, information = information, title = title , description = description, ...)
    lowerfit = gpdfit(-x, -nLower, type = type, information = information, title = title , description = description, ...)
    kernelFit = bkde(x, kernelfit, gridsize=as.integer(N),range.x=c(1.5*min(x),1.5*max(x)))
    if (is.null(title)) title = "GPD Two-Tail Parameter Estimation"
    if (is.null(description)) description = .description()
    new("GPDTAILS",
        call = call,
        method = type,
        kernel = kernelfit,
        data = as.matrix(data),
        threshold = list(upper = nUpper , lower = nLower),
		ptails = list(upper = upper, lower = lower),
		fit = list(upperFit = upperfit, lowerFit = lowerfit, kernelFit = kernelFit),
        title = title,
        description = description)
}