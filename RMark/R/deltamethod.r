#' Compute delta method variance for sum, cumsum, prod and cumprod functions
#' 
#' This function computes the delta method std errors or v-c matrix for a sum,
#' cumsum (vector of cummulative sums), prod (product), or cumprod(vector of
#' cummulative products) of a set of estimates.
#' 
#' This function computes the delta method std errors or v-c matrix for a sum,
#' cumsum (vector of cummulative sums), prod (product), or cumprod(vector of
#' cummulative products).  It uses the function deltamethod from the msm
#' package and constructs the necessary formula for these special cases. It
#' will load the msm pacakge but assumes that it has already been installed.
#' See the msm documentation for a complete description on how the deltamethod
#' function works.  If ses=TRUE, it returns a vector of std errors for each of
#' functions of the estimates contained in mean.  If ses=F, then it returns a
#' v-c matrix for the functions of the estimates contained in mean.  cov is the
#' input v-c matrix of the estimates.
#' 
#' @param function.name Quoted character string of either "sum", "cumsum",
#' "prod" or "cumprod"
#' @param mean vector of estimates used in the function
#' @param cov variance-covariance matrix of the estimates
#' @param ses if TRUE it returns a vector of estimated standard errors of the
#' function of the estimates and if FALSE it returns the variance-covariance
#' matrix
#' @return either a vector of standard errors (ses=TRUE) or a
#' variance-covariance matrix (ses=FALSE)
#' @author Jeff Laake
#' @export
#' @import msm
#' @keywords utility
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' #
#' # The following are examples only to demonstrate selecting different 
#' # model sets for adjusting chat and showing model selection table. 
#' # It is not a realistic analysis.
#' #
#'   data(dipper)
#'   mod=mark(dipper,model.parameters=list(Phi=list(formula=~time)))
#'   rr=get.real(mod,"Phi",se=TRUE,vcv=TRUE)
#'   deltamethod.special("prod",rr$estimates$estimate[1:6],rr$vcv.real)
#'   deltamethod.special("cumprod",rr$estimates$estimate[1:6],rr$vcv.real,ses=FALSE)
#'   deltamethod.special("sum",rr$estimates$estimate[1:6],rr$vcv.real)
#'   deltamethod.special("cumsum",rr$estimates$estimate[1:6],rr$vcv.real,ses=FALSE)
#' }
deltamethod.special=function(function.name,mean,cov,ses=TRUE)
{
#
# This function computes the delta method std errors or v-c matrix
# for a sum, cumsum (vector of cummulative sums), prod (product),
# or cumprod(vector of cummulative products).  It uses the function deltamethod
# from the msm package.  It will load the msm pacakge but assumes that it has already
# been installed.  See the msm documentation for a complete description on how the
# deltamethod function works.  If ses=TRUE, it returns a vector of std errors for each of
# functions of the estimates contained in mean.  If ses=FALSE, then it returns a v-c matrix for the
# functions of the estimates contained in mean.  cov is the input v-c matrix of the estimates.
#
# This function handles the special cases of sum, cumsum, prod, cumprod.  It simply
# constructs the necessary formula or list of formula and passes them onto deltamethod
# and then returns the values.
  if(function.name=="prod")
     return(deltamethod(as.formula(paste("~",paste("x",1:length(mean),sep="",collapse="*"),sep="")),mean,cov,ses))
  if(function.name=="cumprod")
  {
     formula.list=vector("list",length=length(mean))
     formula.list[[1]]=~x1
     for (i in 2:length(mean))
        formula.list[[i]]=as.formula(paste("~",paste("x",1:i,sep="",collapse="*"),sep=""))
     return(deltamethod(formula.list,mean,cov,ses))
  }
  if(function.name=="sum")
     return(deltamethod(as.formula(paste("~",paste("x",1:length(mean),sep="",collapse="+"),sep="")),mean,cov,ses))
  if(function.name=="cumsum")
  {
     formula.list=vector("list",length=length(mean))
     formula.list[[1]]=~x1
     for (i in 2:length(mean))
        formula.list[[i]]=as.formula(paste("~",paste("x",1:i,sep="",collapse="+"),sep=""))
     return(deltamethod(formula.list,mean,cov,ses))
  }
}


