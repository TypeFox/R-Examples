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

# example:
# library(xts)
# data(sp500ret)
# spx = as.xts(sp500ret)
# nn = nrow(spx)
# nx = nn-round(0.9*nn,0)
# if(nx%%50!=0) nx = nx+(50-nx%%50)
# h = (nx/50)-1
# indexin = lapply(1:h, function(j){ seq(1,(nn-nx)+j*50, by=1) })
# indexout = lapply(indexin, function(x){ (tail(x,1)+1):(tail(x,1)+50) })

arfimacv = function(data, indexin, indexout, ar.max = 2, ma.max = 2, criterion = c("rmse","mae","berkowitzp"),
		berkowitz.significance = 0.05, arfima = FALSE, include.mean = NULL, distribution.model = "norm",
		cluster = NULL, external.regressors = NULL, solver = "solnp", solver.control=list(), fit.control=list(),
		return.best=TRUE)
{
	sig = berkowitz.significance
	arnames = paste("ar", 1:ar.max, sep = "")
	manames = paste("ma", 1:ma.max, sep = "")
	.str = NULL
	if(ar.max>0){
		for(i in 1:ar.max){
			.str = c(.str, paste(arnames[i],"=c(0,1),",sep=""))
		}
	}
	if(ma.max>0){
		for(i in 1:ma.max){
			.str = c(.str, paste(manames[i],"=c(0,1),",sep=""))
		}
	}
	if(is.null(include.mean)){
		.str = c(.str, "im = c(0,1)")
	} else{
		.str = c(.str, "im = as.integer(include.mean)")
	}
	if(is.null(arfima)){
		.str = c(.str, ",arf = c(0,1)")
	} else{
		.str = c(.str, ",arf = as.integer(arfima)")
	}
	if(!is.null(external.regressors)){
		.str = c(.str, ",xreg = c(0,1)")
	} else{
		.str = c(.str, ",xreg = 0")
	}
	str = c("d = expand.grid(", paste(.str), ')')
	xstr = paste(str, sep="", collapse="")
	eval(parse(text=xstr))
	# eliminate the zero row
	check = apply(d, 1, "sum")
	if(any(check == 0)){
		idx=which(check==0)
		d = d[-idx,]
	}
	sumar = apply(d[,1:ar.max,], 1, "sum")
	summa = apply(d[,(ar.max+1):(ma.max+ar.max),], 1, "sum")
	# some checks on the indexin and indexout:
	if(!is.list(indexin)) stop("\nindexin must be a list of the training indices points")
	if(!is.list(indexout)) stop("\nindexout must be a list of the testing indices points")
	if(length(indexin)!=length(indexout)) stop("\nlength of indexin list must be the same as length of indexout list")
	if(!is(data, "xts")) stop("\ndata must be an xts matrix")
	if(!is.null(external.regressors) && !is(external.regressors, "xts")) stop("\nexternal.regressors must be an xts matrix")
	if(!is.null(external.regressors) && !all.equal(index(external.regressors), index(data))) stop("\nexternal.regressors time index not equal to data time index")
	m = length(indexin)
	n = dim(d)[1]
	if(!is.null(cluster)){
		clusterEvalQ(cluster, library("rugarch"))
		clusterExport(cluster, c("d", "data", "n", "solver", "external.regressors",  "distribution.model",".zarfimaspec",
						"ar.max", "ma.max", "solver.control", "fit.control","indexin","indexout","sig",".evalcv","nazero"), 
				envir = environment())
		# loop through the indices
		cvlist = vector(mode="list", length=m)
		statsmat = matrix(NA, ncol=3, nrow=nrow(d))
		colnames(statsmat) = c("rmse","mae","berkowitzp")
		evmat = lapply(1:m, function(k) cbind(as.matrix(d), statsmat))
		mxd = nrow(d)
		for(i in 1:m){
			datain = data[indexin[[i]]]
			dataout = data[indexout[[i]]]
			dout = index(dataout)
			din = index(datain)
			dall = unique(c(din,dout))
			clusterExport(cluster,c("dout","din","dall"),envir = environment())
			# loop through the indices
			cvlist[[i]] = parLapply(cluster, as.list(1:n), fun = function(j){
						if(ar.max>0){
							arr = d[j,1:ar.max]
							if(ma.max>0){
								mar = d[j,(ar.max+1):(ma.max+ar.max)]
							} else{
								mar = 0
							}
						} else{
							arr = 0
							if(ma.max>0){
								mar = d[j,1:ma.max]
							} else{
								mar = 0
							}
						}
						if(d[j,"xreg"]==1){
							ex = external.regressors[din]
							exf = external.regressors[dall]
						} else{
							ex = NULL
							exf = NULL
						}
						spec = .zarfimaspec( arOrder = arr, maOrder = mar, 
								include.mean = d[j,'im'], arfima = d[j,'arf'], 
								external.regressors = ex, distribution.model = distribution.model)
						fit = try(arfimafit(spec = spec, data = data[din], solver = solver, 
										solver.control = solver.control, 
										fit.control = fit.control), silent = TRUE)
						while(inherits(fit, 'try-error') || (convergence(fit)==1 | all(is.na(fit@fit$matcoef[,2])))){
							fit = try(arfimafit(spec = spec, data = data[din], solver = "gosolnp", 
											solver.control = solver.control, 
											fit.control = fit.control), silent = TRUE)
						}
						specx = .zarfimaspec( arOrder = arr, maOrder = mar, 
								include.mean = d[j,'im'], arfima = d[j,'arf'], 
								external.regressors = exf, distribution.model = distribution.model)
						setfixed(specx)<-as.list(coef(fit))
						forc = arfimafilter(specx, data = data[dall], n.old=length(din))
						f = fitted(forc)[dout]
						r = data[dout]
						evalx = .evalcv(f, r, coef(fit)["sigma"], skew = nazero(coef(fit)["skew"]), shape=nazero(coef(fit)["shape"]),
								distribution = distribution.model, sig = sig)
						return(evalx)
					})
			evmat[[i]][,c("rmse","mae","berkowitzp")]<-t(sapply(cvlist[[i]], function(x) x))
		}
	} else{
		# loop through the indices
		cvlist = vector(mode="list", length=m)
		statsmat = matrix(NA, ncol=3, nrow=nrow(d))
		colnames(statsmat) = c("rmse","mae","berkowitzp")
		evmat = lapply(1:m, function(k) cbind(as.matrix(d), statsmat))
		mxd = nrow(d)
		for(i in 1:m){
			datain = data[indexin[[i]]]
			dataout = data[indexout[[i]]]
			dout = index(dataout)
			din = index(datain)
			dall = unique(c(din,dout))
			clusterExport(cluster,c("dout","din","dall"),envir = environment())
			# loop through the indices
			for(j in 1:n){
					if(ar.max>0){
						arr = d[j,1:ar.max]
						if(ma.max>0){
							mar = d[j,(ar.max+1):(ma.max+ar.max)]
						} else{
							mar = 0
						}
					} else{
						arr = 0
						if(ma.max>0){
							mar = d[j,1:ma.max]
						} else{
							mar = 0
						}
					}
					if(d[j,"xreg"]==1){
						ex = external.regressors[din]
						exf = external.regressors[dall]
					} else{
						ex = NULL
						exf = NULL
					}
					spec = .zarfimaspec( arOrder = arr, maOrder = mar, 
							include.mean = d[j,'im'], arfima = d[j,'arf'], 
							external.regressors = ex, distribution.model = distribution.model)
					fit = try(arfimafit(spec = spec, data = data[din], solver = solver, 
									solver.control = solver.control, 
									fit.control = fit.control), silent = TRUE)
					while(inherits(fit, 'try-error') || (convergence(fit)==1 | all(is.na(fit@fit$matcoef[,2])))){
						fit = try(arfimafit(spec = spec, data = data[din], solver = "gosolnp", 
										solver.control = solver.control, 
										fit.control = fit.control), silent = TRUE)
					}
					specx = .zarfimaspec( arOrder = arr, maOrder = mar, 
							include.mean = d[j,'im'], arfima = d[j,'arf'], 
							external.regressors = exf, distribution.model = distribution.model)
					setfixed(specx)<-as.list(coef(fit))
					forc = arfimafilter(specx, data = data[dall], n.old=length(din))
					f = fitted(forc)[dout]
					r = data[dout]
					evalx = .evalcv(f, r, coef(fit)["sigma"], skew = nazero(coef(fit)["skew"]), shape=nazero(coef(fit)["shape"]),
							distribution = distribution.model, sig = sig)
					cvlist[[i]] = evalx
			}
			evmat[[i]][,c("rmse","mae","berkowitzp")]<-t(sapply(cvlist[[i]], function(x) x))
		}
	}
	# evaluate across the cross validation sets
	resstats = matrix(NA, ncol=6, nrow=nrow(d))
	colnames(resstats) = c("rmse", "rmse_sd", "mae", "mae_sd", "berkowitzp", "berkowitzp_sd")
	resmat = cbind(as.matrix(d), resstats)
	for(i in 1:nrow(d)){
		resmat[i,"rmse"] = mean(sapply(evmat, function(x) x[i,"rmse"]), na.rm=TRUE)
		resmat[i,"rmse_sd"] = sd(sapply(evmat, function(x) x[i,"rmse"]), na.rm=TRUE)
		resmat[i,"mae"] = mean(sapply(evmat, function(x) x[i,"mae"]), na.rm=TRUE)
		resmat[i,"mae_sd"] = sd(sapply(evmat, function(x) x[i,"mae"]), na.rm=TRUE)
		resmat[i,"berkowitzp"] = mean(sapply(evmat, function(x) x[i,"berkowitzp"]), na.rm=TRUE)
		resmat[i,"berkowitzp_sd"] = mean(sapply(evmat, function(x) x[i,"berkowitzp"]), na.rm=TRUE)
	}
	if(criterion[1]=="berkowitzp"){
		bestmodel = as.numeric(which(resmat[,"berkowitzp"]==max(resmat[,"berkowitzp"])))
		
	} else{
		bestmodel = as.numeric(which(resmat[,criterion[1]]==max(resmat[,criterion[1]])))
	}
	if(return.best){
		# estimate once more:
		if(ar.max>0){
			arr = d[j,1:ar.max]
			if(ma.max>0){
				mar = d[bestmodel,(ar.max+1):(ma.max+ar.max)]
			} else{
				mar = 0
			}
		} else{
			arr = 0
			if(ma.max>0){
				mar = d[bestmodel,1:ma.max]
			} else{
				mar = 0
			}
		}
		if(d[bestmodel,"xreg"]==1){
			ex = external.regressors
		} else{
			ex = NULL
		}
		spec = .zarfimaspec( arOrder = arr, maOrder = mar, 
				include.mean = d[bestmodel,'im'], arfima = d[bestmodel,'arf'], 
				external.regressors = ex, distribution.model = distribution.model)
		fit = try(arfimafit(spec = spec, data = data, solver = solver, solver.control = solver.control, 
						fit.control = fit.control), silent = TRUE)
		while(inherits(fit, 'try-error') | convergence(fit)==1 | all(is.na(fit@fit$matcoef[,2]))){
			fit = try(arfimafit(spec = spec, data = data, solver = "gosolnp", solver.control = solver.control, 
							fit.control = fit.control), silent = TRUE)
		}
	} else{
		bestmodel = NULL
	}
	return(list(bestmodel = fit, cv_matrix = resmat))
}

nazero = function(x)
{
	if(is.na(x)) return(0) else return(x)
}
.evalcv = function(pred, act, sigma, skew=0, shape=0, distribution, sig)
{
	rmse = sqrt(mean( (act-pred)^2 ))
	mae  = mean( abs(act-pred) )
	bk = pdist(distribution, q = act, mu = pred, sigma = sigma, skew = skew, shape = shape)
	nb = qnorm(bk)
	bt = try(BerkowitzTest(data = nb, lags = 1, significance = sig)$LRp, silent=TRUE)
	if(inherits(bt, 'try-error')) bt = NA
	ans = matrix(c(rmse, mae, bt), nrow=1)
	colnames(ans) = c("rmse","mae","berkowitz")
	return(ans)
}