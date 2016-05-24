########### Compute Probability of being under the LOQ, for all observations

#' Compute P(y_ij<LOQ) for all observations
#' 
#' For each observation in the dataset, computes the probability to be below LOQ
#' using the simulations under the model. This is generally performed automatically.
#' 
#' 
#' @usage compute.ploq(npdeObject)
#' @param npdeObject an object returned by a call to \code{\link{npde}} or \code{\link{autonpde}}
#' @return an object of class NpdeObject with a component ploq in the ["results"] slot
#' @author Emmanuelle Comets <emmanuelle.comets@@bichat.inserm.fr>
#' @seealso \code{\link{npde}}, \code{\link{autonpde}}
#' @references K. Brendel, E. Comets, C. Laffont, C. Laveille, and F.
#' Mentre. Metrics for external model evaluation with an application to the
#' population pharmacokinetics of gliclazide. \emph{Pharmaceutical Research},
#' 23:2036--49, 2006.
#' @keywords internal

compute.ploq<-function(npdeObject) {
	# For each observation, computes the probability of being LOQ under the model
	### when the data has different LOQ values, the LOQ is taken to be the smallest LOQ
	### p_LOQ is computed as the % of simulated values under the LOQ
	if(length(npdeObject["data"]["icens"])==0)
		return(npdeObject)
	ploq<-c()
	ploqfull<-rep(NA,dim(npdeObject["data"]["data"])[1])
	tab<-npdeObject["data"]["data"][npdeObject["data"]["not.miss"],] # non-missing data
	tabsim<-npdeObject["sim.data"]["datsim"][npdeObject["data"]["not.miss"],] # corresponding simulated data    
	yobs<-tab[, npdeObject["data"]["name.response"]]
	nrep<-npdeObject["sim.data"]["nrep"]
	ysim<-tabsim$ysim
	idsim<-tabsim$idsim
	
	if(length(npdeObject["data"]["loq"])>0) {
		loq<-npdeObject["data"]["loq"]
		if(npdeObject@options$verbose) cat("Computing p(y<LOQ) using LOQ=",loq,"\n")
	} else {
		yloq<-npdeObject["data"]["data"][npdeObject["data"]["icens"], npdeObject["data"]["name.response"]]
		if(length(unique(yloq))==1) {
			if(npdeObject@options$verbose) cat("Same LOQ for all missing data, loq=",loq,"\n")
			loq<-unique(yloq)
		} else {
			loq<-min(unique(yloq))
			if(npdeObject@options$verbose) cat("Computing p(y<LOQ) for the lowest LOQ, loq=",loq,"\n")
		}
		npdeObject["data"]["loq"]<-loq
	}
	for(isuj in unique(npdeObject["data"]["data"][,"index"])) {
		matsim<-matrix(ysim[idsim==isuj],ncol=nrep)
		tcomp<-apply(matsim,2,"<",loq)
		if(!is.matrix(tcomp)) tcomp<-t(as.matrix(tcomp))
		ycal<-rowMeans(tcomp)
		ploq<-c(ploq,ycal)
	}
	ploqfull[npdeObject["data"]["not.miss"]]<-ploq
	npdeObject["results"]["ploq"]<-ploqfull
	#saving pd
	invisible(npdeObject)
}

########### Compute pd

#' Internal functions used to compute prediction discrepancies
#' 
#' Functions used by \code{npde} and \code{autonpde} to compute prediction
#' discrepancies (not intended for the end-user).
#' 
#' These functions are normally not called by the end-user.
#' 
#' @usage computepd(npdeObject,...)
#' @param npdeObject an object of class NpdeObject
#' @param \ldots additional options to modify
#' @return an object of class NpdeObject; the results slot will now include prediction discrepancies (pd) as well as model predictions (ypred) obtained as the mean of the simulated data for each observation
#' @author Emmanuelle Comets <emmanuelle.comets@@bichat.inserm.fr>
#' @seealso \code{\link{npde}}, \code{\link{autonpde}}
#' @references K. Brendel, E. Comets, C. Laffont, C. Laveille, and F.
#' Mentre. Metrics for external model evaluation with an application to the
#' population pharmacokinetics of gliclazide. \emph{Pharmaceutical Research},
#' 23:2036--49, 2006.
#' @keywords internal

computepd<-function(npdeObject,...) {
	npdeObject@options<-replace.control.options(npdeObject@options,...)
	if(npdeObject["options"]$cens.method=="omit" | length(npdeObject["data"]["icens"])==0)
    object<-computepd.omit(npdeObject) else {
		  if(npdeObject["options"]$cens.method=="cdf") {
		  	if(length(npdeObject["results"]["ploq"])==0)
		  		npdeObject<-compute.ploq(npdeObject)
		    object<-computepd.cdf(npdeObject)
		  }
		  if(npdeObject["options"]$cens.method %in% c("ipred","ppred"))
		    object<-computepd.replace(npdeObject)
		  if(npdeObject["options"]$cens.method=="fix")
		  	object<-computepd.replace(npdeObject,...)
		  if(npdeObject["options"]$cens.method=="loq") {
		  	if(length(npdeObject["data"]["loq"])==0)
		  		npdeObject<-compute.ploq(npdeObject)
		  	object<-computepd.replace(npdeObject)
		  }
    }
  return(object)
}

# Functions to sort the lines/columns of a matrix
linesort<-function(mat) {
  for(i in 1:dim(mat)[1])
    mat[i,]<-sort(mat[i,])
  return(mat)
}
colsort<-function(mat) {
	for(i in 1:dim(mat)[2])
		mat[,i]<-sort(mat[,i])
	return(mat)
}

computepd.omit<-function(npdeObject) {
# Preparing vectors for the different results
	ypredfull<-pdfull<-rep(NaN,dim(npdeObject["data"]["data"])[1])
	ycomp<-npdeObject["data"]["data"][,npdeObject["data"]["name.response"]]
# Extracting only non-missing and non-censored data
	if(length(npdeObject["data"]["icens"])==0) keep<-npdeObject["data"]["not.miss"] else keep<-(npdeObject["data"]["not.miss"] & npdeObject["data"]["data"][,npdeObject["data"]["name.cens"]]==0)
# missing and censored data = NaN in complete data
	ycomp[!keep]<-NaN
	tab<-npdeObject["data"]["data"][keep,] # non-missing data
	tabsim<-npdeObject["sim.data"]["datsim"][keep,] # corresponding simulated data    
	yobs<-tab[,npdeObject["data"]["name.response"]]
	idobs<-tab[,npdeObject["data"]["name.group"]]
	nrep<-npdeObject["sim.data"]["nrep"]
	ysim<-tabsim$ysim
	idsim<-tabsim$idsim
	
	ypredall<-pd<-c()
# Loop on isuj
	for(isuj in unique(idobs)) {
		matsim<-matrix(ysim[idsim==isuj],ncol=nrep)
		ysuj<-yobs[idobs==isuj]
# Compute ypred
		ypred<-rowMeans(matsim)
		ypredall<-c(ypredall,ypred)
# compute pd_ij
		tcomp<-apply(matsim,2,"<",ysuj)
		if(!is.matrix(tcomp)) tcomp<-t(as.matrix(tcomp))
		ycal<-rowMeans(tcomp)
		if(npdeObject["options"]$ties) {
			ycal[!is.na(ycal) & ycal==0]<-1/(2*nrep)
			ycal[!is.na(ycal) & ycal==1]<-1-1/(2*nrep)
		} else {
			idx<-which(!is.na(ycal) & ycal>0 & ycal<1)
			ycal[!is.na(ycal) & ycal==0]<-runif(sum(ycal==0),0,1/nrep)
			ycal[!is.na(ycal) & ycal==1]<-runif(sum(ycal==0),1-1/nrep,1)
			ycal[idx]<-ycal[idx]+runif(length(idx),0,1/nrep)
		}
		pd<-c(pd,ycal)
	}
	pdfull[keep]<-pd
	ypredfull[keep]<-ypredall
	
# Saving results
	if(length(npdeObject["results"]["res"])==0) { # create data.frame res
		res<-data.frame(ypred=ypredfull,ycomp=ycomp,pd=pdfull)
		npdeObject["results"]["res"]<-res
	} else { # append to data.frame res
		npdeObject["results"]["res"]$ypred<-ypredfull
		npdeObject["results"]["res"]$ycomp<-ycomp
		npdeObject["results"]["res"]$pd<-pdfull
	}
	
	return(npdeObject)
}

computepd.replace<-function(npdeObject,fix=NULL) {
# Preparing vectors for the different results
	ypredfull<-pdfull<-rep(NaN,dim(npdeObject["data"]["data"])[1])
	ycomp<-npdeObject["data"]["data"][,npdeObject["data"]["name.response"]]
# Extracting only non-missing data
	keep<-npdeObject["data"]["not.miss"]
# missing data = NaN in complete data
	ycomp[!keep]<-NaN
	tab<-npdeObject["data"]["data"][keep,] # non-missing data
	tabsim<-npdeObject["sim.data"]["datsim"][keep,] # corresponding simulated data    
	yobs<-tab[,npdeObject["data"]["name.response"]]
	icens<-tab[,npdeObject["data"]["name.cens"]]
	idobs<-tab[,npdeObject["data"]["name.group"]]
	if(npdeObject["options"]$cens.method=="fix" & !is.null(fix)) 
		yobs[icens==1]<-fix
	if(npdeObject["options"]$cens.method=="ipred") 
		yobs[icens==1]<-tab[icens==1,npdeObject["data"]["name.ipred"]]
	if(npdeObject["options"]$cens.method=="loq") 
		yobs[icens==1]<-npdeObject["data"]["loq"]
	nrep<-npdeObject["sim.data"]["nrep"]
	ysim<-tabsim$ysim
	idsim<-tabsim$idsim
	
	ypredall<-pd<-c()
# Loop on isuj
	for(isuj in unique(idobs)) {
		matsim<-matrix(ysim[idsim==isuj],ncol=nrep)
		ysuj<-yobs[idobs==isuj]
# Compute ypred
		ypred<-rowMeans(matsim)
		ypredall<-c(ypredall,ypred)
		if(npdeObject["options"]$cens.method=="ppred") ysuj[icens[idobs==isuj]==1]<-ypred[icens[idobs==isuj]==1]
# compute pd_ij
		tcomp<-apply(matsim,2,"<",ysuj)
		if(!is.matrix(tcomp)) tcomp<-t(as.matrix(tcomp))
		ycal<-rowMeans(tcomp)
		if(npdeObject["options"]$ties) {
			ycal[!is.na(ycal) & ycal==0]<-1/(2*nrep)
			ycal[!is.na(ycal) & ycal==1]<-1-1/(2*nrep)
		} else {
			idx<-which(!is.na(ycal) & ycal>0 & ycal<1)
			ycal[!is.na(ycal) & ycal==0]<-runif(sum(ycal==0),0,1/nrep)
			ycal[!is.na(ycal) & ycal==1]<-runif(sum(ycal==0),1-1/nrep,1)
			ycal[idx]<-ycal[idx]+runif(length(idx),0,1/nrep)
		}
		pd<-c(pd,ycal)
	}
	pdfull[keep]<-pd
	ypredfull[keep]<-ypredall
	if(npdeObject["options"]$cens.method=="ppred") 
		ycomp[npdeObject["data"]["icens"]]<-ypredfull[npdeObject["data"]["icens"]]
	if(npdeObject["options"]$cens.method %in% c("ipred","fix","loq")) { 
		ycomp[npdeObject["data"]["icens"]]<-yobs[icens==1]
	}
	# Saving results
	if(length(npdeObject["results"]["res"])==0) { # create data.frame res
		res<-data.frame(ypred=ypredfull,ycomp=ycomp,pd=pdfull)
		npdeObject["results"]["res"]<-res
	} else { # append to data.frame res
		npdeObject["results"]["res"]$ypred<-ypredfull
		npdeObject["results"]["res"]$ycomp<-ycomp
		npdeObject["results"]["res"]$pd<-pdfull
	}
	return(npdeObject)
}

computepd.cdf<-function(npdeObject) {
# Censored observation
### cdf: replace pd_ij with random sample from a uniform distribution [0-p_LOQ_ij], where p_LOQ_ij is the probability of being under the LOQ (==previous computation of pd_ij because y_ij=censoring value)
	
# Preparing vectors for the different results
	ypredfull<-pdfull<-rep(NaN,dim(npdeObject["data"]["data"])[1])
	ycomp<-npdeObject["data"]["data"][,npdeObject["data"]["name.response"]]
# Indicators for censored data to be imputed
	not.miss<-npdeObject["data"]["not.miss"]
# missing data = NaN in complete data
	ycomp[!not.miss]<-NaN
	tab<-npdeObject["data"]["data"][not.miss,] # non-missing data
	tabsim<-npdeObject["sim.data"]["datsim"][not.miss,] # corresponding simulated data    
	yobs<-tab[,npdeObject["data"]["name.response"]]
	ycens<-tab[,npdeObject["data"]["name.cens"]]
	ploq<-npdeObject["results"]["ploq"][not.miss]
	loq<-npdeObject["data"]["loq"]
	has.cens<-(length(loq)>0)
	idobs<-tab[,npdeObject["data"]["name.group"]]
	nrep<-npdeObject["sim.data"]["nrep"]
	ysim<-tabsim$ysim
	idsim<-tabsim$idsim
	ypredall<-pd<-c()
	
# ECO ici comment gerer les indices proprement ?
# Compute prediction discrepancies	
# Loop on isuj
	for(isuj in unique(idobs)) {
		matsim<-matrix(ysim[idsim==isuj],ncol=nrep)
# compute pd_ij
		tcomp<-apply(matsim,2,"<",yobs[idobs==isuj])
	  if(!is.matrix(tcomp)) tcomp<-t(as.matrix(tcomp))
		pdsuj<-ycal<-rowMeans(tcomp)
		ploq.suj<-ploq[idobs==isuj]
# impute pd_ij for censored data
		iobs.loq<-which(ycens[idobs==isuj]==1)
		if(length(iobs.loq)>0)
			pdsuj[iobs.loq]<-runif(length(iobs.loq),0,ploq.suj[iobs.loq])
# Compute ypred
		ypred<-rowMeans(matsim)
  	ypredall<-c(ypredall,ypred)
# If we plan to compute npde afterwards, we need to impute observed & simulated data under the LOQ; for simulated data, we first need to impute pd in the same way
		isim.loq<-which(matsim<loq,arr.ind=TRUE)
		if(npdeObject["options"]$calc.npde & (dim(isim.loq)[1]+length(iobs.loq))>0) {
# when pd=k/K, jitter by imputing an observation to a value between ysim(k) (kth simulated value) and ysim(k+1)
			matsort<-linesort(matsim)
			if(!npdeObject["options"]$ties) {
# Pb of jittering when the imputed pd is 1 (k=K) because no ysim(K+1) 
				ncol<-dim(matsort)[2]
			  matsort<-cbind(matsort,matsort[,ncol]*2-matsort[,(ncol-1)])
			}
			if(length(iobs.loq)>0){
				ids1<-trunc(pdsuj[iobs.loq]*nrep)+1
				if(!npdeObject["options"]$ties) {
					yobs.imp<-yobs[idobs==isuj]
					for(j in 1:length(iobs.loq)) yobs.imp[iobs.loq[j]]<-runif(1, matsort[iobs.loq[j],ids1[j]],matsort[iobs.loq[j],(ids1[j]+1)])
					yobs[idobs==isuj]<-yobs.imp
				} else {
					yobs.imp<-yobs[idobs==isuj]
					yobs.imp[iobs.loq]<-diag(matsort[iobs.loq,ids1,drop=F])
					yobs[idobs==isuj]<-yobs.imp
				}
			}
# saving the imputed data for simulations in the simulated dataset
			if(dim(isim.loq)[1]>0) {
				pdsim.imp<-runif(dim(isim.loq)[1],0,ploq.suj[isim.loq[,1]])
				ids2<-trunc(pdsim.imp*nrep)+1
				isim.repl<-isim.loq
				if(!npdeObject["options"]$ties) {
					ids2[ids2>=nrep]<-nrep-1
					isim.repl[,2]<-ids2
					isim.repl2<-isim.repl
					isim.repl2[,2]<-isim.repl2[,2]+1
					matsim[isim.loq]<-runif(dim(isim.loq)[1], matsort[isim.repl], matsort[isim.repl2])
					} else {
					isim.repl[,2]<-ids2
					matsim[isim.loq]<-matsort[isim.repl]
					}
				ysim[idsim==isuj & ysim<loq]<-matsim[isim.loq]
				}
			}
		pd<-c(pd,pdsuj)
	}

# Dealing with extreme values of pd & smoothing if !(npdeObject["options"]$ties)
	if(npdeObject["options"]$ties) {
		pd[pd==0]<-1/(2*nrep)
		pd[pd==1]<-1-1/(2*nrep)
	} else {
		idx<-which(pd>0 & pd<1)
		pd[pd==0]<-runif(sum(pd==0),0,1/nrep)
		pd[pd==1]<-runif(sum(pd==1),1-1/nrep,1)
		pd[idx]<-pd[idx]+runif(length(idx),0,1/nrep)
	}

	if(npdeObject["options"]$calc.npde) {
		ycomp[not.miss]<-yobs
		npdeObject["sim.data"]["datsim"][not.miss,"ysim.imp"]<-ysim
	}
	ypredfull[not.miss]<-ypredall
	pdfull[not.miss]<-pd
	# Saving results
	if(length(npdeObject["results"]["res"])==0) { # create data.frame res
		res<-data.frame(ypred=ypredfull,ycomp=ycomp,pd=pdfull)
		npdeObject["results"]["res"]<-res
	} else { # append to data.frame res
		npdeObject["results"]["res"]$ypred<-ypredfull
		npdeObject["results"]["res"]$ycomp<-ycomp
		npdeObject["results"]["res"]$pd<-pdfull
	}
    return(npdeObject)
}

########### Compute npde

#' Internal functions used to compute normalised prediction distribution errors (npde)
#' 
#' Functions used by \code{npde} and \code{autonpde} to compute prediction
#' discrepancies (not intended for the end-user).
#' 
#' These functions are normally not called by the end-user.
#' 
#' @usage computenpde(npdeObject,...)
#' @param npdeObject an object of class NpdeObject
#' @param \ldots additional options to modify
#' @return an object of class NpdeObject; the results slot will now include prediction discrepancies (npde) as well as decorrelated observed data, while the sim.data slot will now include decorrelated simulated data
#' @author Emmanuelle Comets <emmanuelle.comets@@bichat.inserm.fr>
#' @seealso \code{\link{npde}}, \code{\link{autonpde}}
#' @references K. Brendel, E. Comets, C. Laffont, C. Laveille, and F.
#' Mentre. Metrics for external model evaluation with an application to the
#' population pharmacokinetics of gliclazide. \emph{Pharmaceutical Research},
#' 23:2036--49, 2006.
#' @keywords internal

computenpde<-function(npdeObject,...) {
	npdeObject@options<-replace.control.options(npdeObject@options,...)
	if(npdeObject["options"]$cens.method=="cdf" & length(npdeObject["data"]["icens"])>0 & length(npdeObject["results"]["res"]$ycomp)==0) {
		cat("With method",npdeObject["options"]$cens.method,"prediction discrepancies need to be computed first \n")
		npdeObject<-computepd(npdeObject)
	}	
	if(npdeObject["options"]$cens.method=="omit" | length(npdeObject["data"]["icens"])==0)
		object<-computenpde.omit(npdeObject) else
			object<-computenpde.loq(npdeObject)
	return(object)
}

computenpde.omit<-function(npdeObject) {
	# Compute (normalised) prediction distribution errors omitting BQL
	# Preparing vectors for the different results
	# Indicators for censored data to be imputed
	if(npdeObject["options"]$cens.method=="omit" & length(npdeObject["data"]["icens"])>0) {
		tab<-npdeObject["data"]["data"]
		not.miss<-(tab[,npdeObject["data"]["name.miss"]]==0 & tab[,npdeObject["data"]["name.cens"]]==0)
	} else not.miss<-npdeObject["data"]["not.miss"]
	tab<-npdeObject["data"]["data"][not.miss,] # non-missing data
	tabsim<-npdeObject["sim.data"]["datsim"][not.miss,] # corresponding simulated data    
	ypred<-pde<-rep(NA,dim(tab)[1])
	yobs.all<-ydobs<-tab[,npdeObject["data"]["name.response"]]
	npdefull<-ydobsfull<-rep(NA,dim(npdeObject["data"]["data"])[1])
	ydsimfull<-rep(NA,dim(npdeObject["sim.data"]["datsim"])[1])
	id<-tab[,npdeObject["data"]["name.group"]]
	nrep<-npdeObject["sim.data"]["nrep"]
	xerr<-0    
	ysim<-ydsim<-tabsim$ysim
	idsim<-tabsim$idsim
	
	# computing pde
	for(isuj in unique(id)) {
		yobs<-yobs.all[id==isuj]
		matsim<-matrix(ydsim[idsim==isuj],ncol=nrep)
		x<-calcnpde(isuj,yobs,matsim,nrep,npdeObject["options"]$decorr.method, npdeObject["options"]$verbose)
		if(x$xerr>0) {
			xerr<-x$xerr
			break
		}
		pde[id==isuj]<-x$pde
		ydsim[idsim==isuj]<-x$ydsim
		ydobs[id==isuj]<-x$ydobs
		ypred[id==isuj]<-x$ypred
	}
	if(xerr>0) {
		cat("The computation of the pde has failed for subject",isuj,"because \n")
		if(xerr==1) {
			if(npdeObject["options"]$decorr.method=="cholesky") cat("the Cholesky decomposition of the covariance matrix of the simulated data could not be obtained.\n")
			if(npdeObject["options"]$decorr.method=="inverse") cat("the covariance matrix of the simulated data could not be diagonalised through eigen().\n")
			if(npdeObject["options"]$decorr.method=="polar") cat("the Cholesky decomposition of the covariance matrix of the simulated data could not be obtained.\n")
		}
		if(xerr==2) cat("the covariance matrix of the simulated data could not be inverted.\n")
		cat("This usually means that the covariance matrix is not positive-definite, or that is is poorly conditioned.\n")
		cat("This can be caused by simulations widely different from observations (in \n")
		cat("other words, a poor model).\n")
		cat("We suggest to plot a prediction interval from the simulated data to check\n")
		cat("whether the simulations are reasonable, and to consider prediction\n")
		cat("discrepancies (obtained without the decorrelation step).\n")
		cat("Prediction discrepancies will now be computed.\n")
		break
	}
	
	# saving pde
	if(xerr==0) {
		if(npdeObject["options"]$ties) {
			pde[pde==0]<-1/(2*nrep)
			pde[pde==1]<-1-1/(2*nrep)
		} else {
			idx<-which(pde>0 & pde<1)
			pde[pde==0]<-runif(sum(pde==0),0,1/nrep)
			pde[pde==1]<-runif(sum(pde==1),1-1/nrep,1)
			pde[idx]<-pde[idx]+runif(length(idx),0,1/nrep)
		}
		npde<-qnorm(pde) 
		npdeObject["results"]["xerr"]<-xerr
		npdefull[not.miss]<-npde
		ydobsfull[not.miss]<-ydobs
		ydsimfull[not.miss]<-ydsim
		if(length(npdeObject["results"]["res"])==0) {
			res<-data.frame(npde=npdefull,ydobs=ydobsfull)
			npdeObject["results"]["res"]<-res
		} else {
			npdeObject["results"]["res"]$ydobs<-ydobsfull
			npdeObject["results"]["res"]$npde<-npdefull
		}
		npdeObject["sim.data"]["datsim"]$ydsim<-ydsimfull
		if(!npdeObject["options"]$calc.pd) { 
			ypredfull<-npdefull
			ypredfull[not.miss]<-ypred
			npdeObject["results"]["res"]$ypred<-ypredfull
		}
	} else {
		npde<-rep(NA,length(pde))
		npdeObject["options"]$calc.pd<-TRUE
	}
	return(npdeObject)
}

computenpde.loq<-function(npdeObject) {
	# Compute (normalised) prediction distribution errors in the presence of BQL data
	
	# ECO TODO: securiser ici, faire test
	if(length(npdeObject["data"]["icens"])>0 & length(npdeObject["results"]["res"]$pd)==0) {
		cat("Please compute prediction discrepancies first \n")
		return(npdeObject)
	}
	# Preparing vectors for the different results
	# Indicators for censored data to be imputed
	not.miss<-npdeObject["data"]["not.miss"]
	tab<-npdeObject["data"]["data"][not.miss,] # non-missing data
	tabsim<-npdeObject["sim.data"]["datsim"][not.miss,] # corresponding simulated data    
	ycens<-tab[,npdeObject["data"]["name.cens"]]
	ypred<-pde<-rep(NA,dim(tab)[1])
	if(npdeObject["options"]$cens.method%in% c("ipred","ppred","cdf")) 		
		yobs.all<-ydobs<-npdeObject["results"]["res"][not.miss,"ycomp"] else 
			yobs.all<-ydobs<-tab[,npdeObject["data"]["name.response"]]
	npdefull<-ydobsfull<-rep(NA,dim(npdeObject["data"]["data"])[1])
	ydsimfull<-rep(NA,dim(npdeObject["sim.data"]["datsim"])[1])
	loq<-npdeObject["data"]["loq"]
	has.cens<-(length(loq)>0)
	id<-tab[,npdeObject["data"]["name.group"]]
	nrep<-npdeObject["sim.data"]["nrep"]
	xerr<-0    
	ysim<-tabsim$ysim
	idsim<-tabsim$idsim
	ypredall<-pd<-c()    
	if(npdeObject["options"]$cens.method=="cdf") ydsim<-tabsim$ysim.imp else ydsim<-ysim
	
	# computing pde
	for(isuj in unique(id)) {
		yobs<-yobs.all[id==isuj]
		matsim<-matrix(ydsim[idsim==isuj],ncol=nrep)
		x<-calcnpde(isuj,yobs,matsim,nrep,npdeObject["options"]$decorr.method, npdeObject["options"]$verbose)
		if(x$xerr>0) {
			xerr<-x$xerr
			break
		}
		pde[id==isuj]<-x$pde
		ydsim[idsim==isuj]<-x$ydsim
		ydobs[id==isuj]<-x$ydobs
		ypred[id==isuj]<-x$ypred
	}
	if(xerr>0) {
		cat("The computation of the pde has failed for subject",isuj,"because \n")
		if(xerr==1) {
			if(npdeObject["options"]$decorr.method=="cholesky") cat("the Cholesky decomposition of the covariance matrix of the simulated data could not be obtained.\n")
			if(npdeObject["options"]$decorr.method=="inverse") cat("the covariance matrix of the simulated data could not be diagonalised through eigen().\n")
			if(npdeObject["options"]$decorr.method=="polar") cat("the Cholesky decomposition of the covariance matrix of the simulated data could not be obtained.\n")
		}
		if(xerr==2) cat("the covariance matrix of the simulated data could not be inverted.\n")
		cat("This usually means that the covariance matrix is not positive-definite, or that is is poorly conditioned.\n")
		cat("This can be caused by simulations widely different from observations (in \n")
		cat("other words, a poor model).\n")
		cat("We suggest to plot a prediction interval from the simulated data to check\n")
		cat("whether the simulations are reasonable, and to consider prediction\n")
		cat("discrepancies.\n")
		cat("Prediction discrepancies will now be computed.\n")
		break
	}
	
	# saving pde
	if(xerr==0) {
		# Dealing with extreme values of npde & smoothing if !(npdeObject["options"]$ties)
		if(npdeObject["options"]$ties) {
			pde[pde==0]<-1/(2*nrep)
			pde[pde==1]<-1-1/(2*nrep)
		} else {
			idx<-which(pde>0 & pde<1)
			pde[pde==0]<-runif(sum(pde==0),0,1/nrep)
			pde[pde==1]<-runif(sum(pde==1),1-1/nrep,1)
			pde[idx]<-pde[idx]+runif(length(idx),0,1/nrep)
		}
		npde<-qnorm(pde) 
		npdeObject["results"]["xerr"]<-xerr
		npdefull[not.miss]<-npde
		npdeObject["results"]["res"]$npde<-npdefull
		ydobsfull[not.miss]<-ydobs
		ydsimfull[not.miss]<-ydsim
		npdeObject["results"]["res"]$ydobs<-ydobsfull
		npdeObject["sim.data"]["datsim"]$ydsim<-ydsimfull
	} else {
		npde<-rep(NA,length(pde))
		npdeObject["options"]$calc.pd<-TRUE
	}
	return(npdeObject)
}


########### Compute npde - decorrelation methods

calcnpde<-function(isuj,yobs,matsim,nrep,decorr.method,verbose) {
	# matsim: simulated data, with BQL data replaced (depending on method)
	# matsim: simulated data, raw
	if (verbose) cat("Computing the npde for subject ",isuj,"\n")
	#Computing decorrelated ysim* and yobs* for subject isuj
	#variance-covariance matrix computed using the cov function
	#Computing ypred
	varsim<-cov(t(matsim))
	moysim<-rowMeans(matsim)
	#computing V-1/2 with Cholesky
	xerr<-0
	if(length(moysim)>1) {
		y<-switch(decorr.method,
							cholesky=decorr.chol(varsim),
							polar=decorr.polar(varsim),
							inverse=decorr.inverse(varsim))
		if(y$xerr==0) ymat<-y$y else xerr<-y$xerr
	} else ymat<-1/sqrt(varsim)
	if(xerr==0) {
		#decorrelation of the simulations
		decsim<-t(ymat)%*%(matsim-moysim)
		decobs<-t(ymat)%*%(yobs-moysim)
		ydsim<-c(decsim)
		#decorrelation of the observations
		ydobs<-decobs
		#Computing the pde
		tcomp<-apply(decsim,2,"<",decobs)
		if(!is.matrix(tcomp)) tcomp<-t(as.matrix(tcomp))
		ycal<-rowMeans(tcomp)
		pde<-ycal
	}
	return(list(xerr=xerr,pde=pde,ydsim=ydsim,ydobs=ydobs,ypred=moysim))
}

decorr.chol<-function(x) {
	xerr<-0
	xmat<-try(chol(x))
	if(is.numeric(xmat)) {
		ymat<-try(solve(xmat))
		if(!is.numeric(ymat)) 
			xerr<-2
	} else 
		xerr<-1
	return(list(y=ymat,xerr=xerr))
}

decorr.inverse<-function(x) {
	xerr<-0
	var.eig<-eigen(x)
	xmat<-try(var.eig$vectors %*% diag(sqrt(var.eig$values)) %*% solve(var.eig$vectors))
	if(is.numeric(xmat)) {
		ymat<-try(solve(xmat))
		if(!is.numeric(ymat)) 
			xerr<-2
	} else 
		xerr<-1
	return(list(y=ymat,xerr=xerr))
}

decorr.polar<-function(x) {
	xerr<-0
	xmat<-try(chol(x))
	if(is.numeric(xmat)) {
		svdec<-svd(xmat)
		umat<-svdec$u %*% t(svdec$v)
		vmat<-t(umat) %*% xmat
		ymat<-try(solve(vmat))
		if(!is.numeric(ymat)) 
			xerr<-2
	} else 
		xerr<-1
	return(list(y=ymat,xerr=xerr))
}

########### Compute distribution of pd/npde under the model, using simulations

#' Compute distribution of pd/npde using simulations
#' 
#' This function is used to built the distribution of pd/npde using the simulations under the model. The default is to build only the distribution of pd, and to sample from N(0,1) when building the distribution of npde under the null hypothesis.
#' 
#' @aliases dist.pred.sim calcnpde.sim
#' @usage dist.pred.sim(npdeObject,nsamp, ...)
#' @param npdeObject an object returned by a call to \code{\link{npde}} or \code{\link{autonpde}}
#' @param nsamp number of datasets (defaults to 100 or to the number of replications if it is smaller)
#' @param \dots additional arguments. Currently only the value of calc.pd and calc.npde may be passed on, and will override their corresponding value in the "options" slot of npdeObject
#' @return an object of class NpdeObject; the ["results"] slot will contain pd and/or npde for a sample of the simulated datasets (depending on whether calc.pd/calc.npde are ), stored in pd.sim and/or npde.sim
#' @author Emmanuelle Comets <emmanuelle.comets@@bichat.inserm.fr>
#' @seealso \code{\link{npde}}, \code{\link{autonpde}}
#' @references K. Brendel, E. Comets, C. Laffont, C. Laveille, and F.
#' Mentre. Metrics for external model evaluation with an application to the
#' population pharmacokinetics of gliclazide. \emph{Pharmaceutical Research},
#' 23:2036--49, 2006.
#' @export
#' @examples
#' 
#' data(theopp)
#' data(simtheopp)
#' x<-autonpde(theopp,simtheopp,1,3,4,boolsave=FALSE)
#' # Use random samples from N(0,1) to obtain a prediction interval on the empirical cdf of the npde
#' plot(x,plot.type="ecdf",bands=TRUE,approx.pi=TRUE)
#' # defaults to computing the pd and npde for 100 simulated datasets (in the theophylline example, this uses all the simulated datasets)
#' x<-dist.pred.sim(x)
#' # Use the npde from the simulated datasets to obtain a prediction interval on the empirical cdf
#' plot(x,plot.type="ecdf",bands=TRUE,approx.pi=FALSE)

dist.pred.sim<-function(npdeObject,nsamp, ...) {
	args1<-match.call(expand.dots=TRUE)
	i1<-match("calc.pd",names(args1))
	calc.pd<-NA
	if(!is.na(i1)) calc.pd<-as.logical(as.character(args1[[i1]]))
	if(is.na(calc.pd)) calc.pd<- npdeObject["options"]$calc.pd
	i1<-match("calc.npde",names(args1))
	calc.npde<-NA
	if(!is.na(i1)) calc.npde<-as.logical(as.character(args1[[i1]])) 
	if(is.na(calc.npde)) calc.npde<- npdeObject["options"]$calc.npde
	if(!calc.pd & !calc.npde) {
		cat("At least one of calc.pd or calc.npde must be TRUE.\n")
		return(npdeObject)		
	}
# ECO not necessary since we're not using the imputed data... why not ??? maybe should ?
# 	if(npdeObject["options"]$cens.method=="cdf" && length(npdeObject["sim.data"]["icens"])>0 && length(npdeObject["sim.data"]["datsim"]$ysim.imp)==0) {
# 		cat("With the cdf method, the imputation step needs to be performed before a call to dist.pred.sim() is made. Please use computepd(x) first.\n")
# 		return(npdeObject)
# 	}
	nrep<-npdeObject["sim.data"]["nrep"]
	if(missing(nsamp)) nsamp<-min(100,nrep)
	# Extracting non-missing data
	keep<-npdeObject["data"]["not.miss"]
	tabsim<-npdeObject["sim.data"]["datsim"][keep,] # simulated data    
#	if(npdeObject["options"]$cens.method=="cdf" && !is.na(grep("ysim.imp",colnames(tabsim)))) yobs<-matrix(tabsim[,"ysim.imp"], ncol=nrep) else yobs<-matrix(tabsim[,"ysim"],ncol=nrep)
# Use original simulations to compute the PI under the model
	idsim<-tabsim$idsim
	ysim<-yobs<-matrix(tabsim[,"ysim"],ncol=nrep)
	colnames(yobs)<-paste("sim",1:nrep,sep="")
	if(nsamp<nrep) {
		isample<-sort(sample(1:nrep,nsamp))
		yobs<-yobs[,isample]
	}
	idobs<-npdeObject["data"]["data"][keep,npdeObject["data"]["name.group"]]
	
	pd<-npde<-matrix(nrow=0,ncol=nsamp,dimnames=list(NULL,colnames(yobs)))
	# Loop on isuj
	for(isuj in unique(idobs)) {
		matsim<-matrix(ysim[idsim==isuj],ncol=nrep)
		ysuj<-yobs[idobs==isuj,]
		if(calc.pd) {
			pdsuj<-c()
			# compute pdsim_ij
			for(i in 1:nsamp) {
				tcomp<-apply(matsim,2,"<",ysuj[,i])
				if(!is.matrix(tcomp)) tcomp<-t(as.matrix(tcomp))
				ycal<-rowMeans(tcomp)
# pdsim_ij will be within the sequence seq(0,1-1/nrep,1/nrep) by construction				
				pdsuj<-c(pdsuj,ycal)
			}
# If ties: add 1/(2*nrep) to center the distribution to 1/2*nrep ; 1-1/2*nrep
# else add U(0,1/nrep)
			if(npdeObject["options"]$ties)
				pdsuj<-pdsuj+1/(2*nrep)
			else
				pdsuj<-pdsuj+runif(length(pdsuj),0,1/nrep)
			pd<-rbind(pd,matrix(pdsuj,ncol=nsamp))
		}
		if(calc.npde) {
			y<-calcnpde.sim(ysuj,matsim,nrep,npdeObject["options"]$decorr.method)
			if(y$xerr==0) {
				pde<-y$pde
				if(npdeObject["options"]$ties)
					pde<-pde+1/(2*nrep)
				else
					pde<-pde+runif(length(pde),0,1/nrep)
				npde<-rbind(npde,qnorm(pde))
			} else {
				cat("Problem computing npde for subject",isuj,"\n")
				return(npdeObject)
			}
		}
	}
	
	# Saving results
	if(calc.pd) npdeObject["results"]["pd.sim"]<-pd
	if(calc.npde) npdeObject["results"]["npde.sim"]<-npde	
	return(npdeObject)	
}

calcnpde.sim<-function(yobs,matsim,nrep,decorr.method) {
	# yobs: matrix with data treated as observed (1 column=1 set of data=1 'subject')
	# matsim: simulated data to be used to decorrelate
	# returns decorrelated yobs* for each 'subject'
	#variance-covariance matrix computed depending on "method"
	varsim<-cov(t(matsim))
	moysim<-rowMeans(matsim)
	#computing V-1/2
	xerr<-0
	if(length(moysim)>1) {
		y<-switch(decorr.method,
							cholesky=decorr.chol(varsim),
							polar=decorr.polar(varsim),
							inverse=decorr.inverse(varsim))
		if(y$xerr==0) ymat<-y$y else xerr<-y$xerr
	} else ymat<-1/sqrt(varsim)
	pde<-yobs
	if(xerr==0) {
		#decorrelation of the simulations & the observations
		decsim<-t(ymat)%*%(matsim-moysim)
		decobs<-t(ymat)%*%(yobs-moysim)
		#Computing the pde
		for(i in 1:dim(yobs)[2]) {
			tcomp<-apply(decsim,2,"<",decobs[,i])
			if(!is.matrix(tcomp)) tcomp<-t(as.matrix(tcomp))
			ycal<-rowMeans(tcomp)
			pde[,i]<-ycal
		}
	}
	return(list(xerr=xerr,pde=pde))
}

####################################################################################
