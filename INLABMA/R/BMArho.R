#require(multicore)

#
#Compute marginal for rho using marginal likelihhod x prior
#Compute the marginals for all the other parameters using BMA


#
#Version of splinefun to return 0 outside x-range
#
mysplinefun<-function(x, y = NULL,
               method = c("fmm", "periodic", "natural", "monoH.FC")[1],
               ties = mean)
{
	xmin<-min(x)
	xmax<-max(x)

	ff<-splinefun(x=x, y=y, method=method, ties=ties)

	fff<-function(x, deriv)
	{

		sapply(x, function(X){
		if(X<xmin |X>xmax)
			return(0)
		else
			return(ff(X))
		})
	}

}



#Compute and re-scale marginal using splines
#x: x-valuies 
#logy: log(y)-values (posterior LOG-density)
#logp: log-prior for x-values
#usenormal: Use Normal approximation??
#

fitmarg<-function(x, logy, logp=0, usenormal=FALSE)
{
	if(!usenormal)
	{
	#func = splinefun(x, logy-max(logy))
	#post = exp(logp+func(x))

	logpost<-logy-max(logy)+logp-max(logp)
	post<-exp(logpost)

	post.func = splinefun(x, post)

	## normalize
	xx = seq(min(x), max(x),  len = 1000)
	#z = sum(post.func(xx)) * diff(xx)[1]
	z=integrate(post.func, min(x), max(x))$value
	post.func = mysplinefun(x=xx, y=(post.func(xx) / z) )

	}
	else
	{
			
	xx = seq(min(x), max(x),  len = 1000)
	meany=sum(x*exp(logp+logy))/sum(exp(logp+logy))
	sdy=sqrt(   sum(((x-meany)^2)*exp(logp+logy))/sum(exp(logp+logy))  )

	post.func = function(x){dnorm(x, mean=meany, sd=sdy)}

	}

	return(post.func)
} 


#
#Fits a marginal using BMA
#margs: list of matrices with the marginals are in INLA (i.e., 2-col matrix)
#ws: weights

fitmargBMA<-function(margs, ws, len=100)
{

	ws<-ws/sum(ws)

	xmin<-max(unlist(lapply(margs, function(X){min(X[,1])})))
	xmax<-min(unlist(lapply(margs, function(X){max(X[,1])})))


	xx<-seq(xmin, xmax, len=len)


	margsws<-lapply(1:length(margs), function(i){
		func<-fitmarg(margs[[i]][,1], margs[[i]][,2])
		ws[i]*func(xx)
	})

	margsws<-do.call(cbind, margsws)

	d<-data.frame(x=xx, y=apply(margsws, 1, sum))
	names(d)<-c("x", "y")

	return(d)
}


#
#Weighted sum of summary results in matrices
#
fitmatrixBMA<-function(models, ws, item)
{
	lmatrix<-lapply(models, function(X){X[[item]]})

	auxbma<-ws[1]*lmatrix[[1]]
        for(i in 2:length(lmatrix)){auxbma<-auxbma+ws[i]*lmatrix[[i]]}

	return(auxbma)
}

#
#Weighted sum of summary results in lists
#
fitlistBMA<-function(models, ws, item)
{
	nlist<-names(models[[1]][[item]])
	auxlist<-as.list(rep(NA, length(nlist)))
	names(auxlist)<-nlist

	for(ele in nlist)
	{
		lmatrix<-lapply(models, function(X){X[[item]][[ele]]})

		auxbma<-ws[1]*lmatrix[[1]]
	        for(i in 2:length(lmatrix)){auxbma<-auxbma+ws[i]*lmatrix[[i]]}

		auxlist[[ele]]<-auxbma
	}
	return(auxlist)
}



#
#Weighted sum of marginals 
#
#The elements to handle are lists with different marginals (in matrix form)
#
fitmargBMA2<-function(models, ws, item)
{

	if(is.null(models[[1]][[item]]))
		return(NULL)

	nlist<-names(models[[1]][[item]])
	auxlist<-as.list(rep(NA, length(nlist)))
	names(auxlist)<-nlist

	for(ele in nlist)
	{
		lmatrix<-lapply(models, function(X){X[[item]][[ele]]})
		#xxr<-range(unlist(lapply(lmatrix, function(X){X[,1]})))

		#Compute range using a weighted sum of min and max
		xxr<-c(NA, NA)
		xxr[1]<-sum(ws*unlist(lapply(lmatrix,function(X){min(X[,1])})))
		xxr[2]<-sum(ws*unlist(lapply(lmatrix,function(X){max(X[,1])})))

		xx<-seq(xxr[1], xxr[2], length.out=81)

		auxbma<-rep(0, length(xx))
	        for(i in 1:length(lmatrix)){
			auxspl<- mysplinefun(x=lmatrix[[i]][,1], y=lmatrix[[i]][,2])
			auxbma<-auxbma+ws[i]*auxspl(xx)
			}

		auxlist[[ele]]<-cbind(xx, auxbma)
	}
	return(auxlist)
}


#Returns a summary of the fitted.values using BMA
#
#models: list of INLA models
#rrho: Vector of values of rho (same as number of models)
#logrhoprior:vector of logdensities of prior for rho

BMArho<-function(models, rho, logrhoprior=rep(1, length(rho)) )
{
	mlik<-unlist(lapply(models, function(X){X$mlik[1]}))

	post.func<-fitmarg(rho, mlik, logrhoprior)

#	mlik.func = splinefun(rho, mlik - max(mlik))
#
#	post = exp(logrhoprior) * exp(mlik.func(rho))
#	post.func = splinefun(rho, post)
#
#	## normalize
#	rrho = seq(min(rho), max(rho),  len = 1000)
#	z = sum(post.func(rrho)) * diff(rrho)[1]
#	post.func = splinefun(rrho, post.func(rrho) / z)

	#Weights for BMA
	ws<-(post.func(rho))
	ws<-ws/sum(ws)

	#Fitted values with BMA

	#Compute BMA
	fvalues<-mclapply(1:length(models), function(X){ws[X]*models[[X]]$summary.fitted.values[,1]})
	fvalues<-data.frame(fvalues)
	fvalues<-apply(fvalues, 1, sum)

	return(fvalues)
}

#Elements worth BMA
#
# "summary.fixed" "marginals.fixed"
# "summary.lincomb" "marginals.lincomb"
# "summary.random" "marginals.random"
# "summary.linear.predictor" "marginals.linear.predictor"
# "summary.fitted.values" "marginals.fitted.values"

#Performs BMA on a number of things
#
#
#
INLABMA<-function(models, rho, logrhoprior=rep(1, length(rho)), impacts=FALSE, usenormal=FALSE )
{

	#require(INLA)

	mlik<-unlist(lapply(models, function(X){X$mlik[1]}))
	post.func<-fitmarg(rho, mlik, logrhoprior, usenormal)

	#Weights for BMA
	ws<-(post.func(rho))
	ws<-ws/sum(ws)

	mfit<-list(rho=list())
	mfit$rho$marginal<-data.frame(x=seq(min(rho), max(rho), len=100))
	mfit$rho$marginal$y<-post.func(mfit$rho$marginal$x)
	mfit$rho$marginal<-as.matrix(mfit$rho$marginal)

#	mfit$rho$mean<-sum(mfit$rho$marginal$x*mfit$rho$marginal$y)*diff(mfit$rho$marginal$x)[1]
#	mfit$rho$sd<-sqrt(sum(((mfit$rho$marginal$x-mfit$rho$mean))^2*mfit$rho$marginal$y)*diff(mfit$rho$marginal$x)[1])

	#Compute some quantiles
	#fquant<-function(x, ff, qtile){integrate(ff, 0, x)$value-qtile}
	#qvect<-c(.025, .5, .975)
	#mfit$rho$quantiles<-sapply(qvect, 
	#   function(qtile){
	#	uniroot(fquant, interval=c(min(rho),max(rho)), ff=post.func, qtile=qtile)$root
	#   }
	#)
	#names(mfit$rho$quantiles)<-as.character(qvect)

	#Summary statistics
	margsum <- INLA::inla.zmarginal(mfit$rho$marginal, TRUE)
	mfit$rho$mean<-margsum$mean
	mfit$rho$sd<-margsum$sd
	mfit$rho$quantiles<-unlist(margsum[-c(1:2)])



	#Results which are stored in matrices
	mateff<-c("summary.fixed", "summary.lincomb", 
	#"summary.random", 
		"summary.linear.predictor", "summary.fitted.values")
	#summary.fixed
	#mfit$summary.fixed<-fitmatrixBMA(lapply(models, function(X){X$summary.fixed}), ws)

	lmat<-mclapply(mateff, function(X){fitmatrixBMA(models, ws, X)})
	names(lmat)<-mateff

	mfit<-c(mfit, lmat)

	#Results stored in lists
	listeff<-c("dic", "cpo")

	leff<-mclapply(listeff, function(X){fitlistBMA(models, ws, X)})
	names(leff)<-listeff

	mfit<-c(mfit, leff)

	#Results are marginals
	listmarg<-c("marginals.fixed", "marginals.lincomb",
		"marginals.lincomb.derived", #"marginals.random",
		"marginals.linear.predictor", #"marginals.fitted.values",
		"marginals.hyperpar", "marginals.spde2.blc")

	margeff<-mclapply(listmarg, function(X){fitmargBMA2(models, ws, X)})
	names(margeff)<-listmarg

	mfit<-c(mfit, margeff)


	#Impacts
	mfit$impacts<-FALSE
	if(impacts)
	{
		mfit$impacts<-TRUE

		summimp<-c("summary.total.impacts", "summary.direct.impacts", 
				"summary.indirect.impacts")
		matsummimp<-mclapply(summimp, function(X){fitmatrixBMA(models, ws, X)})
		names(matsummimp)<-summimp
		mfit<-c(mfit, matsummimp)

		margimp<-c("marginals.total.impacts","marginals.direct.impacts",
				"marginals.indirect.impacts")
		lmargimp<-mclapply(margimp, function(X){fitmargBMA2(models, ws, X)})
		names(lmargimp)<-margimp

		mfit<-c(mfit, lmargimp)
	
		#Recompute impacts summaries
		mfit<-recompute.impacts(mfit)

	}

	return(mfit)
}
