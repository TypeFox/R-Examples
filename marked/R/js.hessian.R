#' Compute variance-covariance matrix for fitted JS model
#' 
#' A wrapper function that sets up call to hessian function to compute and then
#' invert the hessian.
#' 
#' @param model fitted JS model from function crm
#' @export
#' @return variance-covariance matrix for specified model or the model
#' object with the stored vcv depending on whether the model has already been run
#' @author Jeff Laake <jeff.laake@@noaa.gov>
js.hessian=function(model)
{
	object=NULL
#   Previously run model object
	if(!is.null(model$results))
	{
		object=model
		model=model$results
		if(model$options$accumulate)
		{
			capture.output(model_data<-js.accumulate(object$data$data,model$model_data,
							object$data$nocc,object$data$freq,chunk_size=model$options$chunk_size))
		} else
			model_data=model$model_data
		scale=set.scale(names(object$model.parameters),model_data,1)
	} else
#   Called within model fitting code
	{
		scale=model$scale
		model_data=model$model_data
	}	
#nobstot number of unique caught at least once by group if applicable
	markedfunc_eval=0
	jsenv=environment()
	vcv=numDeriv::hessian(js.lnl,scale.par(model$beta,scale),model_data=model$model_data,nobstot=model$ns,jsenv=jsenv)
	vcv=try(solvecov(vcv))
	if(class(vcv)[1]=="try-error")
	{
		warning("Unable to invert hessian")
		return(NULL)
	}
	scale=unlist(scale)
	vcv=vcv$inv/outer(scale,scale,"*")
	colnames(vcv)=names(scale)
	rownames(vcv)=names(scale)
	if(is.null(object))
		return(vcv)
	else
	{
		model$beta.vcv=vcv
		object$results=model
		return(object)
	}
	
}

