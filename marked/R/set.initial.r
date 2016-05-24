#' Set initial values 
#' 
#' Sets initial values specified in a list.
#' 
#' @param pars character vector of parameter names
#' @param dml list of design matrices created by \code{\link{create.dm}} from
#' formula and design data
#' @param initial list of vectors for parameter initial values
#' @return List of initial values for each parameter in the model
#' @author Jeff Laake <jeff.laake@@noaa.gov>
set.initial=function(pars,dml,initial)
{
#   if this is a previously run model, get initial values from it
	if(class(initial)[1]=="crm")
		if(class(initial)[2]=="mcmc")
			initial=lapply(initial$results$beta,function(x){z=x$mean;names(z)=rownames(x);z})	    
	    else
			initial=initial$results$beta
#   create par vector from initial filling in 0 for initial value if not all
#   values or any values are specified
	par=vector("list",length(pars))
	names(par)=pars
#	par=c(par,initial)
	ptype=NULL
#   if no initial list, create one with NULL values
	if(is.null(initial))
	{
		initial=vector("list",length=length(pars))
		names(initial)=pars
	}
	for(parx in pars)
	{
		init=initial[[parx]]
		if(is.null(dml[[parx]]$fe))
		{
			par[[parx]]=NULL
			next
		}
		if(is.null(init))
		{
			par[[parx]]=rep(0,ncol(dml[[parx]]$fe))
		} else
		{
			if(length(init)==1 &is.null(names(init)))
				par[[parx]]=c(init,rep(0,ncol(dml[[parx]]$fe)-1))
			else
			{
				if(is.null(names(init)))
				{
					if(length(init)!=ncol(dml[[parx]]$fe))
						stop(paste("For",parx,",length of initial vector does not match number of parameters."))
					else
						par[[parx]]=init
				} else
				{
					beta.names=colnames(dml[[parx]]$fe)
					par[[parx]]=rep(0,length(beta.names))
					par[[parx]][beta.names%in%names(init)]=init[which(names(init)%in%beta.names)]
				}
			}
		}
		ptype=c(ptype,rep(parx,ncol(dml[[parx]]$fe)))	
		names(par[[parx]])=colnames(dml[[parx]]$fe)
	}
	return(list(par=par,ptype=ptype))
}
