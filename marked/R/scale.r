#' Scaling functions
#' 
#' Set scale, scale dm and scale/unscale parameters 
#'
#' @usage 	set.scale(pars,model_data,scale)
#'  
#'	        scale.dm(model_data,scale)
#' 
#'          scale.par(par,scale)
#'
#'          unscale.par(par,scale)
#' 
#' @aliases set.scale scale.dm scale.par unscale.par
#' @param pars character vector of parameter names
#' @param par list of parameter vectors or vector of parameter values
#' @param scale list or vector of parameter scales
#' @param model_data list of data/design objects
#' @return List of scale values for set.scale, model.data with scaled design matrices for scale.dm,
#' vector of scaled parameter values for scale.par, and list of unscaled parameter vectors for unscale.par
#' @author Jeff Laake <jeff.laake@@noaa.gov>
set.scale=function(pars,model_data,scale)
{
	scale.list=vector("list",length(pars))
	names(scale.list)=pars
	if(!is.null(scale)&&!is.list(scale)&&all(scale==1))
	{
		for(parx in pars)
			scale.list[[parx]]=rep(1,ncol(model_data[[paste(parx,".dm",sep="")]]))
	}else
	{
		for(parx in pars)
		{
			if(is.null(scale[[parx]]))
				scale.list[[parx]]=apply(model_data[[paste(parx,".dm",sep="")]],2,function(x) mean(x[x!=0]))
			else
			{
				if(length(scale[[parx]])==1)
				   scale.list[[parx]]=rep(scale[[parx]],ncol(model_data[[paste(parx,".dm",sep="")]]))
			    else
				   if(length(scale[[parx]])!=ncol(model_data[[paste(parx,".dm",sep="")]]))
					   stop(paste("For",parx,"length of scale does not match length of parameters\n"))
				   else
					   scale.list[[parx]]=scale[[parx]]				   
			}		
		}
	}
	for(parx in pars)
	   names(scale.list[[parx]])=colnames(model_data[[paste(parx,".dm",sep="")]])	
	return(scale.list)
}
scale.dm=function(model_data,scale)
{
	pars=names(scale)
	for(parx in pars)
		model_data[[paste(parx,".dm",sep="")]]=t(t(as.matrix(model_data[[paste(parx,".dm",sep="")]]))/scale[[parx]])
    return(model_data)
}
scale.par=function(par,scale)
{
	pars=names(scale)
	for(parx in pars)
		par[[parx]]=par[[parx]]*scale[[parx]]
    return(unlist(par,use.names=FALSE))
}
unscale.par=function(par,scale)
{
	pars=names(scale)
	snames=factor(unlist(sapply(names(scale),function(x) rep(x,length(scale[[x]])))),levels=pars)
	par.list=split(par,snames)
	for(parx in pars)
	{
		names(par.list[[parx]])=names(scale[[parx]])
		par.list[[parx]]=par.list[[parx]]/scale[[parx]]
	}
    return(par.list)
}
