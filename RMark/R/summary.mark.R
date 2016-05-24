#' Summary of MARK model parameters and results
#' 
#' Creates a summary object of either a MARK model input or model output which
#' includes number of parameters, deviance, AICc, the beta and real parameter
#' estimates and optionally standard errors, confidence intervals and
#' variance-covariance matrices.  If there are several groups in the data, the
#' output is structured by group.
#' 
#' The structure of the summary of the real parameters depends on the type of
#' model and the value of the argument \code{se} and \code{showall}.  If
#' \code{se=F} then only the estimates of the real parameters are shown and
#' they are summarized the result element \code{reals} in PIM format.  The
#' structure of \code{reals} depends on whether the PIMS are upper triangular
#' ("Triang") or a row ("Square" although not really square). For the upper
#' triangular format, the values are passed back as a list of matrices where
#' the list is a list of parameter types (eg Phi and p) and within each type is
#' a list for each group containing the pim as an upper triangular matrix
#' containing the real parameter estimate.  For square matrices, \code{reals}
#' is a list of matrices with a list element for each parameter type, but there
#' is not a second list layer for groups because in the returned matrix each
#' group is a row in the matrix of real estimates.  If \code{se=TRUE} then
#' estimates, standard error (se), lower and upper confidence limits (lcl, ucl)
#' and a "Fixed" indicator is passed for each real parameter.  If the pims for
#' the model were simplified to represent the unique real parameters (unique
#' rows in the design matrix), then it is possible to restict the summary to
#' only the unique parameters with \code{showall=FALSE}.  This argument only
#' has an affect if \code{se=TRUE}. If \code{showall=FALSE}, \code{reals} is
#' returned as a dataframe of the unique real parameters specified in the
#' model.  This does not mean they will all have unique values and it includes
#' all "Fixed" real parameters and any real parameters that cannot be
#' simplified in the case of parameters such as "pent" in POPAN or "Psi" in
#' "Multistrata" that use the multinomial logit link. Use of
#' \code{showall=FALSE} is of limited use but provided for completeness.  In
#' most cases the default of \code{showall=TRUE} will be satisfactory.  In this
#' case, \code{reals} is a list of dataframes with a list element for each
#' parameter type.  The dataframe contains the estimate, se,lcl, ucl,fixed and
#' the associated default design data for that parameter (eg time,age, cohort
#' etc).  The advantage of retrieving the reals in this format is that it is
#' the same regardless of the model, so it enables model averaging the real
#' parameters over different models with differing numbers of unique real
#' parameters.
#' 
#' @usage \method{summary}{mark}(object,...,se=FALSE,vc=FALSE,showall=TRUE,show.fixed=FALSE,brief=FALSE)
#' @param object a MARK model object
#' @param se if FALSE the real parameter estimates are output in PIM format
#' (eg. triangular format); if TRUE, they are displayed as a list with se and
#' confidence interval
#' @param vc if TRUE the v-c matrix of the betas is included
#' @param showall if FALSE it only returns the values of each unique parameter
#' value
#' @param show.fixed if FALSE, each fixed value given NA; otherwise the fixed
#' real value is used. If se=TRUE, default for show.fixed=TRUE
#' @param brief if TRUE, does not show real parameter estimates
#' @param ... additional non-specified argument for S3 generic function
#' @export 
#' @return A list with each of the summarized objects that depends on the
#' argument values. Only the first 4 are given if it is a summary of a model
#' that has not been run. \item{model}{type of model (e.g., CJS)}
#' \item{title}{user define title if any} \item{model.name}{descriptive name of
#' fitted model} \item{call}{call to make.mark.model used to construct the
#' model} \item{npar}{number of fitted parameters} \item{lnl}{-2xLog Likelihood
#' value} \item{npar}{Number of parameters (always the number of columns in
#' design matrix)} \item{chat}{Value of over-dispersion constant if not equal
#' to 1} \item{npar.unadjusted}{number of estimated parameters from MARK if
#' different than npar } \item{AICc}{Small sample corrected AIC using npar;
#' named qAICc if chat not equal to 1} \item{AICc.unadjusted}{Small sample
#' corrected AIC using npar.unadjusted; prefix of q if chat not equal to 1}
#' \item{beta}{dataframe of beta parameters with estimate, se, lcl, ucl}
#' \item{vcv}{variance-covariance matrix for beta} \item{reals}{list of lists,
#' dataframes or matrices depending on value of se and the type of model
#' (triangular versus square PIMS) (see details above)}
#' @author Jeff Laake
#' @export summary.mark print.summary.mark coef.mark
#' @keywords utility
summary.mark <-function(object,...,se=FALSE,vc=FALSE,showall=TRUE,show.fixed=FALSE,
                    brief=FALSE)
{
# -------------------------------------------------------------------------------------------------------------
#
# summary.mark  - creates a summary of either a MARK model input or model output
#
# Arguments:
#   object    - a MARK model object
#   se        - if FALSE the real parameter estimates are output in PIM (eg. triangular format); if TRUE, they
#                are displayed as a list with se and confidence interval
#   vc        - if TRUE the v-c matrix of the betas is printed
#   showall   - if FALSE only unique real parameters are sumamrized; added 10 Jan 06
#   show.fixed- if TRUE fixed values are returned rather than NA in place of fixed values
#                 if se=TRUE, show.fixed default is TRUE, otherwise FALSE
#   brief     - if TRUE, does not show reals
#
# Value:
#   a list of summary values
#
# Functions used: get.real, setup.parameters
#
# -------------------------------------------------------------------------------------------------------------
#
#
  model=load.model(object)
  if(se & missing(show.fixed))show.fixed=TRUE
#
# Display baseline info about model (type, name, title etc)
#
x=list(model=model$model,title=model$title,model.name=model$model.name,model.call=model$call)
#
# If displaying model input only show call 
#
if(is.null(model$output))
{
  class(x)="summary.mark"
  return(x)
}
#
# If summarizing model output, show num of parameters, deviance, AICc
#
x$npar=model$results$npar
if(!is.null(model$results$npar.unadjusted))
   x$npar.unadjusted=model$results$npar.unadjusted
   x$lnl=model$results$lnl
   x$AICc=model$results$AICc
if(!is.null(model$results$npar.unadjusted))
   x$AICc.unadjusted=model$results$AICc.unadjusted
if(is.null(model$chat))
    chat=1
else
    chat=model$chat
if(chat!=1)
{
   K=model$results$npar
   qaicc= model$results$lnl/chat + 2*K + 2*K*(K+1)/(model$results$n-K-1)
   x$chat=chat
   x$qAICc=qaicc
   if(!is.null(model$results$AICc.unadjusted))
    {
       K=model$results$npar.unadjusted
       qaicc.unadjusted= model$results$lnl/chat + 2*K + 2*K*(K+1)/(model$results$n-K-1)
       x$qAICc.unadjusted=qaicc.unadjusted
    }
}
#
# Display beta coefficients and optionally its v-c matrix
#
  x$beta=coef(model)
if(vc)
{
    vcv=model$results$beta.vcv*chat
    row.names(vcv)=row.names(model$results$beta)
    colnames(vcv)=row.names(vcv)
    x$vcv=vcv
 }
if(!brief)
{
#
# For each parameter type in the model, display the real parameters (by group if any) as either a list
# or in PIM format (se=FALSE)
#
   parameters=model$parameters
   parameter.names=names(parameters)
   x$reals=vector("list",length=length(parameter.names))
   for(i in 1:length(parameter.names))
   {
      x$reals[[i]]=get.real(model,parameter.names[i],se=se,show.fixed=show.fixed)
      if(se)
      {
        if(!showall)
             x$reals[[i]]=x$reals[[i]][!duplicated(x$reals[[i]]$par.index),,drop=FALSE]
        if(!show.fixed)x$reals[[i]]=x$reals[[i]][x$reals[[i]]$fixed!="Fixed",]
      }
      if(parameters[[i]]$type%in%c("Triang","STriang") && parameters[[i]]$pim.type%in%c("time","age")&&!se)
        for (j in 1:length(x$reals[[i]]))
        {
           row.matrix=matrix(x$reals[[i]][[j]]$pim[1,],nrow=1)
           colnames(row.matrix)=colnames(x$reals[[i]][[j]]$pim)
           rownames(row.matrix)=model$begin.time[min(j,length(model$begin.time))]
           x$reals[[i]][[j]]$pim=row.matrix
        }
   }
   names(x$reals)=parameter.names
}
x$brief=brief
class(x)="summary.mark"
return(x)
}
#' MARK model beta parameters
#' 
#' @usage \method{coef}{mark}(object,...)
#' @param object a MARK model object
#' @param ... additional non-specified argument for S3 generic function
#' @export coef.mark
#' @export
#' @return A vector or dataframe of beta estimates
#' @author Jeff Laake
#' @export
#' @keywords utility
coef.mark=function(object,...)
{
   model=load.model(object)
   if(is.null(model$chat))
       chat=1
   else
       chat=model$chat
   if(chat==1)
       beta=model$results$beta
   else
   {
       beta=model$results$beta
       beta$se=sqrt(chat)*beta$se
       beta$lcl=beta$estimate - 1.96*beta$se
       beta$ucl=beta$estimate + 1.96*beta$se
   }
   return(beta)
}


