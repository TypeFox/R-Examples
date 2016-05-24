#' Compute HMM matrices
#' 
#' Computes gamma, dmat and delta (initial) matrices(arrays) and returns them in a list.
#' 
#' @param object fitted crm model (must be an HMM model)
#' @param ddl design data list
#' @return list with gamma, dmat and delta arrays
#' @author Jeff Laake
#' @export compute_matrices
#' @keywords utility
compute_matrices=function(object,ddl=NULL)
{
	if(!substr(object$model,1,3)=="HMM"&!substr(object$model,1,4)=="MVMS")
	{
		message("Not an HMM model. Returning NULL")
		return(NULL)
	}
	if(is.null(ddl))
		ddl=make.design.data(object$data,object$design.parameters)
	ddl=set.fixed(ddl,object$model.parameters)
	if(substr(object$model,1,4)=="MVMS")
	{
		obslevels=object$data$ObsLevels
		sup=object$data$fct_sup(list(obslevels=obslevels))
	} else
		sup=NULL
	if(is.null(object$data$strata.list)|substr(object$model,1,4)=="MVMS"){
		mx=object$data$m
	}else{
		mx=list(ns=length(object$data$strata.list$states),na=length(object$data$strata.list[[names(object$data$strata.list)[names(object$data$strata.list)!="states"]]]))
	}
	dml=create.dml(ddl,model.parameters=object$model.parameters,design.parameters=object$design.parameters,chunk_size=object$results$options$chunk_size)
	mat=HMMLikelihood(par=object$results$beta,xx=object$data$ehmat,mx=mx,T=object$data$nocc,xstart=object$data$start,freq=object$data$freq,
			fct_dmat=object$data$fct_dmat,fct_gamma=object$data$fct_gamma,fct_delta=object$data$fct_delta,ddl=ddl,dml=dml,parameters=object$model.parameters,sup=sup,
			return.mat=TRUE)
	if(!is.null(object$data$strata.labels))
		state.names=c(object$data$strata.labels,"Dead")
	else
		state.names=c("Alive","Dead")
	obs.names=object$data$ObsLevels
	dimnames(mat$gamma)[3:4]=c(list(state.names),list(state.names))
	names(dimnames(mat$gamma))=c("Id","Occasion","From_state","To_state")
	dimnames(mat$dmat)[3:4]=c(list(obs.names),list(state.names))
	names(dimnames(mat$dmat))=c("Id","Occasion","Observation","State")
	return(mat)
}
#' Computes backward probabilities 
#' 
#' Computes backward probability sequence for a set of capture histories
#' 
#' @param object fitted crm model (must be an HMM model)
#' @param ddl design data list
#' @author Jeff Laake
#' @return array of backward probabilities (one for each id, state, occasion)
#' @export backward_prob
#' @keywords utility
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. See page 61.
#' @examples
#' #
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # cormack-jolly-seber model
#' data(dipper)
#' mod=crm(dipper,model="hmmcjs")
#' backward_prob(mod)
#' }

backward_prob=function(object,ddl=NULL)
{  	
	if(!substr(object$model,1,3)=="HMM"&!substr(object$model,1,4)=="MVMS")
	{
		message("Not an HMM model. Returning NULL")
		return(NULL)
	}
	if(!is.null(object$results$mat))
	{
		dmat=object$results$mat$dmat
		gamma=object$results$mat$gamma
	}else
	{
		matlist=compute_matrices(object=object,ddl=ddl)
		dmat=matlist$dmat
		gamma=matlist$gamma
	}
	x=object$data$ehmat
	T=object$data$nocc
	first=object$data$start[,2]
	m=object$data$m
	beta=array(NA,dim=c(nrow(x),ncol(x),m))
	# Loop over capture histories
	for(i in 1:nrow(x))
	{
		occ=T
		beta[i,occ,]=rep(1,m)
		# Loop over occasions for this encounter history (x)
		for(t in T:(first[i]+1))
		{
			occ=occ-1
			# Compute backward probability for this occasion
			beta[i,occ,]=gamma[i,t-1,,]%*%diag(dmat[i,t,x[i,t],])%*%beta[i,occ+1,]  
		}
	}
	return(beta)
} 
#' Local decoding of HMM 
#' 
#' Computes state predictions one at a time for each occasion for each individual
#' 
#' @param object fitted crm model (must be an HMM model)
#' @param ddl design data list
#' @param state.names names for states used to label output; if NULL uses strata.labels + Dead state
#' @author Jeff Laake
#' @return matrix of state predictions
#' @export local_decode
#' @keywords utility
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. See page 80.
#' @examples
#' #
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # cormack-jolly-seber model
#' data(dipper)
#' mod=crm(dipper,model="hmmcjs")
#' local_decode(mod)
#' }

local_decode=function(object,ddl=NULL,state.names=NULL)
{  	
	if(!substr(object$model,1,3)=="HMM"&!substr(object$model,1,4)=="MVMS")
	{
		message("Not an HMM model. Returning NULL")
		return(NULL)
	}
	if(is.null(state.names))
		if(!is.null(object$data$strata.labels))
			state.names=c(object$data$strata.labels,"Dead")
		else
			state.names=c("Alive","Dead")
	result=loglikelihood(object,ddl=ddl)
	result$beta=backward_prob(object,ddl=ddl)
	stateprob=result$alpha*result$beta/exp(result$lnl)
	states=apply(stateprob,c(1,2),function(x){ if(any(is.na(x))) return(NA) else return(state.names[which.max(x)])})
	return(states)
} 
#' Global decoding of HMM 
#' 
#' Computes sequence of state predictions for each individual
#' 
#' @param object fitted crm model (must be an HMM model)
#' @param ddl design data list; will be computed if NULL
#' @param state.names names for states used to label output; if NULL uses strata.labels + Dead state
#' @author Jeff Laake
#' @return matrix of state predictions
#' @export global_decode
#' @keywords utility
#' @references Zucchini, W. and I.L. MacDonald. 2009. Hidden Markov Models for Time Series: An Introduction using R. Chapman and Hall, Boca Raton, FL. See page 82.
#' @examples
#' #
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # cormack-jolly-seber model
#' data(dipper)
#' mod=crm(dipper,model="hmmcjs")
#' global_decode(mod)
#' }
global_decode=function(object,ddl=NULL,state.names=NULL)
{  	
	if(!substr(object$model,1,3)=="HMM"&!substr(object$model,1,4)=="MVMS")
	{
		message("Not an HMM model. Returning NULL")
		return(NULL)
	}
	if(is.null(state.names))
		if(!is.null(object$data$strata.labels))
			state.names=c(object$data$strata.labels,"Dead")
		else
			state.names=c("Alive","Dead")
	parmlist=compute_matrices(object,ddl=ddl)
	dmat=parmlist$dmat
	gamma=parmlist$gamma
	delta=parmlist$delta
	x=object$data$ehmat
	T=object$data$nocc
	m=object$data$m
	first=object$data$start[,2]
	states=matrix(NA,nrow=nrow(x),ncol=T)
	for(i in 1:nrow(x))
	{
		psi=matrix(NA,nrow=T,ncol=m)
		# Assign psi value at first occasion
		psi[first[i],]=delta[i,]%*%diag(dmat[i,first[i],x[i,first[i]],]) 
		for(t in (first[i]+1):T)
		{
			psi[t,]=apply(psi[t-1,]*gamma[i,t-1,,],2,max)%*%diag(dmat[i,t,x[i,t],])
		}
		states[i,T]=which.max(psi[T,])
		for(t in (T-1):first[i])
			states[i,t]=which.max(psi[t,]*gamma[i,t,,states[i,t+1]])
	}
	states=t(apply(states,1,function(x) state.names[x]))
	return(states)
}
