#' Print model results 
#' 
#' Provides a printed simple summary of the model results.
#' 
#' @usage \method{print}{crm}(x,...)
#' @param x crm model result or list of model results
#' @param ... generic arguments not used here
#' @return prints a simple summary of the model to the screen and returns NULL. 
#' @author Jeff Laake
#' @seealso \code{\link{crm}}
#' @keywords utility
#' @export
#' @method print crm

print.crm=function(x,...)
{
   if(mode(x)=="character")x=load.model(x)
   if(!is.null(x$results))x=x$results
   if(class(x)[2]=="admb" & class(x)[3]=="cjs")
   {
	   class(x)[1]="admb"
       print(x)
   }
   else
   {
	   cat("\ncrm Model Summary\n")
       if(class(x)[2]=="mcmc")
	       cat("\nNpar : ",sum(sapply(x$beta,nrow)))   
       else
       {
	       cat("\nNpar : ",sum(sapply(x$beta,length)))
	       cat("\n-2lnL: ",x$neg2lnl)
	       cat("\nAIC  : ",x$AIC)
       }
       cat("\n\nBeta\n")
       print(coef(x))
   }
   invisible(x)
}
#' Print model table from model list
#' 
#' @usage \method{print}{crmlist}(x,...)
#' @param x list of model results
#' @param ... generic arguments not used here
#' @return None
#' @author Jeff Laake
#' @export
#' @method print crmlist
#' @seealso \code{\link{crm}}
#' @keywords utility
print.crmlist<-function(x,...)
{
	if(!is.null(x$model.table))
		print(x$model.table)
	else 
		cat("No model.table is available")
	invisible()
}
