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


fgarchKappa<-function(lambda, delta, eta1, eta2, fk, ghlambda, shape, skew, cond.density,...)
{
	kappa = try(expr=integrate(.ffunE, lower = -Inf, upper = Inf, lambda, delta, eta1, eta2, fk, ghlambda, shape, skew, cond.density,...)[[1]],
			silent=TRUE)
	if(inherits(kappa, "try-error")){
		kappa<-NA}	
	kappa
}

.ffunE<-function(x, lambda, delta, eta1, eta2, fk, ghlambda, shape, skew, cond.density,...)
{   # A function implemented by Alexios Ghalanos
	# Compute Expectation Value
	kdelta=delta+fk*lambda
	cond.density = cond.density[1]
	if (cond.density == "norm"){
		fun = (((abs(x - eta2) - eta1*(x - eta2)))^kdelta) * dnorm(x)
	}
	else if(cond.density == "ged") {
		fun =  (((abs(x - eta2) - eta1*(x - eta2)))^kdelta) * dged(x, shape = shape)
	}
	else if(cond.density == "std") {
		fun =  (((abs(x - eta2) - eta1*(x - eta2)))^kdelta) * dstd(x, shape = shape)
	}
	else if(cond.density == "snorm") {
		fun =  (((abs(x - eta2) - eta1*(x - eta2)))^kdelta) * dsnorm(x, skew = skew)
	}
	else if(cond.density == "sged") {
		fun =  (((abs(x - eta2) - eta1*(x - eta2)))^kdelta) * dsged(x, shape = shape, skew = skew)
	}
	else if(cond.density == "sstd") {
		fun =  (((abs(x - eta2) - eta1*(x - eta2)))^kdelta) * dsstd(x, shape = shape, skew = skew)
	}
	else if(cond.density == "nig") {
		fun =  (((abs(x - eta2) - eta1*(x - eta2)))^kdelta) * dsnig(x, shape = shape, skew = skew)
	}
	else if(cond.density == "ghyp") {
		fun =  (((abs(x - eta2) - eta1*(x - eta2)))^kdelta) * dsgh(x, shape = shape, skew = skew, lambda = ghlambda)
	}
	else if(cond.density == "jsu") {
		fun =  (((abs(x - eta2) - eta1*(x - eta2)))^kdelta) * djsu(x, skew = skew, shape = shape)
	}
	else if(cond.density == "ghst") {
		fun =  (((abs(x - eta2) - eta1*(x - eta2)))^kdelta) * dsghst(x, skew = skew, shape = shape)
	}
	else{
		temp<-paste("d",cond.density,sep="")
		.ddist<-eval(parse(text=paste(temp)))
		fun = (((abs(x - eta2) - eta1*(x - eta2)))^kdelta) * .ddist(x,...)
	}
	# Return Value:
	fun
}

gjrgarchKappa<-function(gm, ghlambda, shape, skew, cond.density,...)
{
	kappa = try(expr=integrate(.gjrfunE, lower = -Inf, upper = Inf, gm, ghlambda, shape, skew, cond.density,...)[[1]],
			silent=TRUE)
	if(inherits(kappa, "try-error")){
		kappa<-NA}
	else{kappa<-integrate(.gjrfunE, lower = -Inf, upper = Inf, gm, ghlambda, shape, skew, cond.density,...)[[1]]}
	
	kappa
}

.gjrfunE<-function(x, gm, ghlambda, shape, skew, cond.density,...)
{   # A function implemented by Alexios Ghalanos
	# Compute Expectation Value
	cond.density = cond.density[1]
	if (cond.density == "norm"){
		fun = (x^2 + gm*(x^2)*(x<0)) * dnorm(x)
	}
	else if(cond.density == "ged") {
		fun =  (x^2 + gm*(x^2)*(x<0)) * dged(x, shape = shape)
	}
	else if(cond.density == "std") {
		fun =  (x^2 + gm*(x^2)*(x<0)) * dstd(x, shape = shape)
	}
	else if(cond.density == "snorm") {
		fun =  (x^2 + gm*(x^2)*(x<0)) * dsnorm(x, skew = skew)
	}
	else if(cond.density == "sged") {
		fun =  (x^2 + gm*(x^2)*(x<0)) * dsged(x, shape = shape, skew = skew)
	}
	else if(cond.density == "sstd") {
		fun =  (x^2 + gm*(x^2)*(x<0)) * dsstd(x, shape = shape, skew = skew)
	}
	else if(cond.density == "nig") {
		fun =  (x^2 + gm*(x^2)*(x<0)) * dsnig(x, shape = shape, skew = skew)
	}
	else if(cond.density == "ghyp") {
		fun =  (x^2 + gm*(x^2)*(x<0)) * dsgh(x, shape = shape, skew = skew, lambda = ghlambda)
	}
	else if(cond.density == "jsu") {
		fun =  (x^2 + gm*(x^2)*(x<0)) * djsu(x, skew = skew, shape = shape)
	}
	else if(cond.density == "ghst") {
		fun =  (x^2 + gm*(x^2)*(x<0)) * dsghst(x, skew = skew, shape = shape)
	}
	else{
		temp<-paste("d",cond.density,sep="")
		.ddist<-eval(parse(text=paste(temp)))
		fun =  (x^2 + gm*(x^2)*(x<0)) * .ddist(x,...)
	}
	# Return Value:
	fun
}


# probability that x<0
pneg<-function(ghlambda, shape, skew, cond.density,...)
{
	kappa = try(expr=integrate(.pnegfunE, lower = -Inf, upper = 0, ghlambda, shape, skew, cond.density,...)[[1]],
			silent=TRUE)
	if(inherits(kappa, "try-error")){
		kappa<-NA}
	else{kappa<-integrate(.pnegfunE, lower = -Inf, upper = 0, ghlambda, shape, skew, cond.density,...)[[1]]}
	
	kappa
}

.pnegfunE<-function(x, ghlambda, shape, skew, cond.density,...)
{   # A function implemented by Alexios Ghalanos
	# Compute Expectation Value
	cond.density = cond.density[1]
	if (cond.density == "norm"){
		fun = dnorm(x)
	}
	else if(cond.density == "ged") {
		fun = dged(x, shape = shape)
	}
	else if(cond.density == "std") {
		fun = dstd(x, shape = shape)
	}
	else if(cond.density == "snorm") {
		fun =dsnorm(x, skew = skew)
	}
	else if(cond.density == "sged") {
		fun = dsged(x, shape = shape, skew = skew)
	}
	else if(cond.density == "sstd") {
		fun = dsstd(x, shape = shape, skew = skew)
	}
	else if(cond.density == "nig") {
		fun = dsnig(x, shape = shape, skew = skew)
	}
	else if(cond.density == "ghyp") {
		fun = dsgh(x, shape = shape, skew = skew, lambda = ghlambda)
	}
	else if(cond.density == "jsu") {
		fun = djsu(x, skew = skew, shape = shape)
	}
	else if(cond.density == "ghst") {
		fun = dsghst(x, skew = skew, shape = shape)
	}
	else{
		temp<-paste("d",cond.density,sep="")
		.ddist<-eval(parse(text=paste(temp)))
		fun = .ddist(x,...)
	}
	# Return Value:
	fun
}

egarchKappa<-function(ghlambda, shape, skew, cond.density,...)
{
	kappa = try(expr=integrate(.efunE, lower = -Inf, upper = Inf, ghlambda, shape, skew, cond.density,...)[[1]],
			silent=TRUE)
	if(inherits(kappa, "try-error")){
		kappa<-NA}
	else{kappa<-integrate(.efunE, lower = -Inf, upper = Inf, ghlambda, shape, skew, cond.density,...)[[1]]}
	
	kappa
}

.efunE<-function(x, ghlambda, shape, skew, cond.density,...)
{   # A function implemented by Alexios Ghalanos
	# Compute Expectation Value
	cond.density = cond.density[1]
	if(cond.density == "norm"){
		fun = abs(x) * dnorm(x)
	}
	else if(cond.density == "ged") {
		fun = abs(x) * dged(x, shape = shape)
	}
	else if(cond.density == "std") {
		fun = abs(x) * dstd(x, shape = shape)
	}
	else if(cond.density == "snorm") {
		fun = abs(x) * dsnorm(x, skew = skew)
	}
	else if(cond.density == "sged") {
		fun = abs(x) * dsged(x, shape = shape, skew = skew)
	}
	else if(cond.density == "sstd") {
		fun = abs(x) * dsstd(x, shape = shape, skew = skew)
	}
	else if(cond.density == "nig") {
		fun = abs(x) * dsnig(x, shape = shape, skew = skew)
	}
	else if(cond.density == "ghyp") {
		fun = abs(x) * dsgh(x, shape = shape, skew = skew, lambda = ghlambda)
	}
	else if(cond.density == "jsu") {
		fun = abs(x) * djsu(x, skew = skew, shape = shape)
	}
	else if(cond.density == "ghst") {
		fun = abs(x) * dsghst(x, skew = skew, shape = shape)
	}
	else{
		temp<-paste("d",cond.density,sep="")
		.ddist<-eval(parse(text=paste(temp)))
		fun = abs(x) * .ddist(x,...)
	}
	# Return Value:
	fun
}

aparchKappa<-function(gm, delta, ghlambda, shape, skew, cond.density,...)
{
	kappa = try(expr=integrate(.afunE, lower = -Inf, upper = Inf, gm, delta, ghlambda, shape, skew, cond.density,...)[[1]],
			silent=TRUE)
	if(inherits(kappa, "try-error")){
		kappa<-NA}
	else{kappa<-integrate(.afunE, lower = -Inf, upper = Inf, gm, delta, ghlambda, shape, skew, cond.density,...)[[1]]}
	
	kappa
}

.afunE<-function(x, gm, delta, ghlambda, shape, skew, cond.density,...)
{   # A function implemented by Alexios Ghalanos
	# Compute Expectation Value
	cond.density = cond.density[1]
	if (cond.density == "norm"){
		fun = ((abs(x)-gm*x)^(delta)) * dnorm(x)
	}
	else if(cond.density == "ged") {
		fun =  ((abs(x)-gm*x)^(delta)) * dged(x, shape = shape)
	}
	else if(cond.density == "std") {
		fun =  ((abs(x)-gm*x)^(delta)) * dstd(x, shape = shape)
	}
	else if(cond.density == "snorm") {
		fun =  ((abs(x)-gm*x)^(delta)) * dsnorm(x, skew = skew)
	}
	else if(cond.density == "sged") {
		fun =  ((abs(x)-gm*x)^(delta)) * dsged(x, shape = shape, skew = skew)
	}
	else if(cond.density == "sstd") {
		fun =  ((abs(x)-gm*x)^(delta)) * dsstd(x, shape = shape, skew = skew)
	}
	else if(cond.density == "nig") {
		fun =  ((abs(x)-gm*x)^(delta)) * dsnig(x, shape = shape, skew = skew)
	}
	else if(cond.density == "ghyp") {
		fun =  ((abs(x)-gm*x)^(delta)) * dsgh(x, shape = shape, skew = skew, lambda = ghlambda)
	}
	else if(cond.density == "jsu") {
		fun =  ((abs(x)-gm*x)^(delta)) * djsu(x, skew = skew, shape = shape)
	}
	else if(cond.density == "ghst") {
		fun =  ((abs(x)-gm*x)^(delta)) * dsghst(x, skew = skew, shape = shape)
	}
	else{
		temp<-paste("d",cond.density,sep="")
		.ddist<-eval(parse(text=paste(temp)))
		fun =  ((abs(x)-gm*x)^(delta)) * .ddist(x,...)
	}
	# Return Value:
	fun
}