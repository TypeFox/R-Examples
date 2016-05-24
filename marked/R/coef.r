#' Extract coefficients
#' 
#' Extracts the beta coefficients from the model results.
#' 
#' @usage \method{coef}{crm}(object,...)
#' @param object crm model result
#' @param ... generic arguments not used here
#' @return  returns a dataframe with estimates and standard
#' errors and confidence intervals if hessian=TRUE on model run.
#' @author Jeff Laake
#' @export 
#' @seealso \code{\link{crm}}
#' @keywords utility
coef.crm=function(object,...)
{
	if("results"%in%names(object)) object=object$results
	if(class(object)[2]=="mcmc")
	{
		beta=do.call(rbind,object$beta)
		indices=grep("\\.",rownames(beta))
		rownames(beta)[-indices]=paste(rownames(beta)[-indices],"(Intercept)",sep=".")
	}
	else
	{
		if(class(object)[2]=="admb" & class(object)[3]=="cjs")
		{
			class(object)[1]="admb"
			beta=coef(object)
		}
		else
		{
			beta=data.frame(Estimate=unlist(object$beta))
			if(!is.null(object$beta.vcv))
			{
				beta$se=sqrt(diag(object$beta.vcv[1:length(beta$Estimate),1:length(beta$Estimate)]))
				beta$lcl=beta$Estimate - 1.96*beta$se
				beta$ucl=beta$Estimate + 1.96*beta$se
			}
		}
	}
	return(beta)
}
