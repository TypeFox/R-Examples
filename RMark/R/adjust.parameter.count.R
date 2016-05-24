#' Adjust count of estimated parameters
#' 
#' Modifies number of estimated parameters and the resulting AICc value for
#' model selection.
#' 
#' When a model is run the parameter count determined by MARK is stored in
#' \code{results$npar} and the AICc value is stored in \code{results$AICc}.  If
#' the argument \code{adjust} is set to TRUE in the call to
#' \code{\link{run.mark.model}} and MARK determined that the design matrix was
#' not full rank (i.e., the parameter count is less than the columns of the
#' design matrix), then the parameter count from MARK is stored in
#' \code{results$npar.unadjusted} and AICc in \code{results$AICc.unadjusted}
#' and \code{results$npar} is set to the number of columns of the design matrix
#' and \code{results$AICc} uses the assumed full rank value of \code{npar}.
#' This function allows the parameter count to be reset to any value less than
#' or equal to the number of columns in the design matrix.  If
#' \code{results$npar.unadjusted} exists it is kept as is.  If it doesn't
#' exist, then the current values of \code{results$npar} and
#' \code{results$AICc} are stored in the \code{.unadjusted} fields to maintain
#' the values from MARK, and the new adjusted values defined by the function
#' argument \code{npar} are stored in \code{results$npar} and
#' \code{results$AICc}. In the example below, the CJS model Phi(t)p(t) is
#' fitted with the call to \code{\link{mark}} which defaults to
#' \code{adjust=TRUE}.  This is used to show how \code{adjust.parameter.count}
#' can be used to adjust the count to 11 from the full rank count of 12.
#' Alternatively, the argument \code{adjust=FALSE} can be added to prevent the
#' adjustment which is appropriate in this case because Phi(6) and p(6) are
#' confounded.
#' 
#' @param model MARK model object
#' @param npar Value of count of estimated parameters
#' @return model: the mark model object with the adjustments made
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{run.mark.model}},\code{\link{model.table}}
#' @keywords utility
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(dipper)
#' ptime=list(formula=~time)
#' Phitime=list(formula=~time)
#' dipper.phitime.ptime=mark(dipper,model.parameters=list(Phi=Phitime, p=ptime))
#' dipper.phitime.ptime=adjust.parameter.count(dipper.phitime.ptime,11)
#' dipper.phitime.ptime=mark(dipper,model.parameters=list(Phi=Phitime, p=ptime),
#'                            adjust=FALSE)
#' }
adjust.parameter.count <-
function(model,npar)
{
#
# adjust.parameter.count - modifies count of estimated parameters and resulting AICc calculation
#
# Arguments:
#
#  model - mark model object
#  npar  - number of estimated parameters 
#
# Value:
#
#  model - mark model object
#
model=load.model(model)
if(model$results$npar!=npar)
{
   if(npar<=dim(model$results$beta)[1])
   {
      message("\nNumber of parameters adjusted from ",model$results$npar," to ", npar,"\n") 
      if(is.null(model$results$AICc.unadjusted))
      {
         model$results$AICc.unadjusted=model$results$AICc
         model$results$npar.unadjusted=model$results$npar
         AICc.unadjusted=model$results$AICc
      }    
      else
         AICc.unadjusted=model$results$AICc.unadjusted
      model$results$npar=npar
      AICc=model$results$lnl+ 2*npar +2*npar*(npar+1)/(model$results$n - npar -1)
      model$results$AICc=AICc
      message("Adjusted AICc = ",AICc,"\nUnadjusted AICc = ", AICc.unadjusted,"\n") 
  }
  else
  {
     message("Specified number of parameters > number of beta parameters. No adjustment made.\n")
  }
}
else
  message("Specified number of parameters matches current value. No adjustment made.\n")
return(model)}
