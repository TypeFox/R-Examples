#' Prints summary of MARK model parameters and results
#'  
#' @usage   \method{print}{summary.mark}(x,...)
#' @param x list resulting from call to \code{summary}
#' @param ... additional non-specified argument for S3 generic function
#' @export
#' @return None
#' @author Jeff Laake
#' @export print.summary.mark
#' @seealso summary.mark
#' @keywords utility
print.summary.mark <-
function(x,...)
{
#
# Display baseline info about model (type, name, title etc)
#
if(is.null(x$npar))
   stype="Input"
else
   stype="Output"
cat(paste(stype,"summary for",x$model, "model\n"))
if(x$title!="")cat("Title:",x$title,"\n")
cat("Name :",x$model.name,"\n")
#
# If displaying model input only show call 
#
if(stype=="Input")
{
  cat("Call : ")
  print(x$model.call)
  cat("\n")
}
else
#
# IF displaying model output, show num of parameters, deviance, AICc
#
{
 if(is.null(x$npar.unadjusted))
   cat("\nNpar : ",x$npar)
 else
   cat("\nNpar : ",x$npar,paste(" (unadjusted=",x$npar.unadjusted,")",sep=""))
 cat("\n-2lnL: ",x$lnl)
 if(is.null(x$AICc.unadjusted))
    cat("\nAICc : ",x$AICc)
 else
   cat("\nAICc : ",x$AICc,paste(" (unadjusted=",x$AICc.unadjusted,")",sep=""))
 if(!is.null(x$chat))
 {
    cat("\nchat : ",x$chat)
    if(is.null(x$AICc.unadjusted))
       cat("\nQAICc: ",x$qAICc)
    else
       cat("\nQAICc: ",x$qAICc,paste(" (unadjusted=",x$qAICc.unadjusted,")",sep=""))
 }
#
# Display beta coefficients and optionally its v-c matrix
#
 cat("\n\nBeta\n")
 print(x$beta)
 if(!is.null(x$vcv))
 {
    cat("\nV-C matrix for beta\n")
    print(x$vcv)
 }
if(x$brief)return()
#
# For each parameter type in the model, display the real parameters (by group if any) as either a list
# or in PIM format (se=FALSE)
#
 if(is.data.frame(x$reals))
 {
    cat("\n\nReal Parameters\n")
    print(x$reals)
 }
 else
 {
    parameter.names=names(x$reals)
    for(i in 1:length(parameter.names))
    {
       cat("\n\nReal Parameter",parameter.names[i])
       if(is.data.frame(x$reals[[i]]))
       {
           cat("\n")
           print(x$reals[[i]][,1:7])
       }else
       {
          if(is.list(x$reals[[i]]))
             for(j in 1:length(x$reals[[i]]))
             {
                cat("\n")
                cat(names(x$reals[[i]])[j],"\n")
                if(is.list(x$reals[[i]][[j]]))
                   print(x$reals[[i]][[j]]$pim,na.print="")
                else
                   print(x$reals[[i]][[j]])
             }
          else
          {
             cat("\n")
             print(x$reals[[i]])
          }
       }
    }
 }
}
invisible()
}
