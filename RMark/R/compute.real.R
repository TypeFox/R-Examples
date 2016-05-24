#' Compute estimates of real parameters
#' 
#' Computes real estimates and var-cov from design matrix (design) and
#' coefficients (beta) using specified link functions
#' 
#' The estimated real parameters can be derived from the estimated beta
#' parameters, a completed design matrix, and the link function specifications.
#' MARK produces estimates of the real parameters, se and confidence intervals
#' but there are at least 2 situations in which it is useful to be able to
#' compute them after running the analysis in MARK: 1) adjusting confidence
#' intervals for estimated over-dispersion, and 2) making estimates for
#' specific values of covariates.  The first case is done in
#' \code{\link{get.real}} with a call to this function.  It is done by
#' adjusting the estimated standard error of the beta parameters by multiplying
#' it by the square root of \code{chat} to adjust for over-dispersion.  A
#' normal 95% confidence interval is computed for the link estimate (estimate
#' +/- 1.96*se) and this is then back-transformed to the real parameters using
#' \code{\link{inverse.link}} with the appropriate inverse link function for
#' the parameter to construct a 95% confidence interval for the real parameter.
#' There is one exception. For parameters using the \code{mlogit}
#' transformation, a \code{logit} transformation of each individual real Psi
#' and its se are used to derive the confidence interval. The estimated
#' standard error for the real parameter is also scaled by the square root of
#' the over-dispersion constant \code{chat} stored in \code{model$chat}. But,
#' the code actually computes the variance-covariance matrix rather than
#' relying on the values from the MARK output because real estimates will
#' depend on any individual covariate values used in the model which is the
#' second reason for this function.
#' 
#' New values of the real parameter estimates can easily be computed by simply
#' changing the values of the covariate values in the design matrix and
#' computing the inverse-link function using the beta parameter estimates.  The
#' covariate values to be used can be specified in one of 2 ways. 1) Prior to
#' making a call to this function, use the functions
#' \code{\link{find.covariates}} to extract the rows of the design matrix with
#' covariate values and either fill in those values aautomatically with the
#' options provided by \code{\link{find.covariates}} or edit those values to be
#' the ones you want and then use \code{\link{fill.covariates}} to replace the
#' values into the design matrix and use it as the value for the argument
#' \code{design}, or 2) automate this step by specifying a value for the
#' argument \code{data} which is used to take averages of the covariate values
#' to fill in the covariate entries of the design matrix.  In computing real
#' parameter estimates from individual covariate values it is important to
#' consider the scale of the individual covariates. By default, an analysis
#' with MARK will standardize covariates by subtracting the mean and dividing
#' by the standard deviation of the covariate value. However, in the
#' \code{RMark} library all calls to MARK.EXE do not standardize the covariates
#' and request real parameter estimates based on the mean covariate values.
#' This was done because there are many instances in which it is not wise to
#' use the standardization implemented in MARK and it is easy to perform any
#' standardization of the covariates with R commands prior to fitting the
#' models.  Also, with pre-standardized covariates there is no confusion in
#' specifying covariate values for computation of real estimates.  If the model
#' contains covariates and the argument \code{design} is not specified, the
#' design matrix is extracted from \code{model} and all individual covariate
#' values are assigned their mean value to be consistent with the default in
#' the MARK analysis.
#' 
#' If a value for \code{beta} is given, those values are used in place of the
#' estimates \code{model$results$beta$estimate}.
#' 
#' @param model MARK model object
#' @param beta estimates of beta parameters for real parameter computation
#' @param design design matrix for MARK model
#' @param data dataframe with covariate values that are averaged for estimates
#' @param se if TRUE returns std errors and confidence interval of real
#' estimates
#' @param vcv logical; if TRUE, sets se=TRUE and returns v-c matrix of real
#' estimates
#' @return A data frame (\code{real}) is returned if \code{vcv=FALSE};
#' otherwise, a list is returned also containing vcv.real: \item{real}{ data
#' frame containing estimates, and if se=TRUE or vcv=TRUE it also contains
#' standard errors and confidence intervals and notation of whether parameters
#' are fixed or at a boundary} \item{vcv.real}{variance-covariance matrix of
#' real estimates}
#' @author Jeff Laake
#' @export
#' @seealso
#' \code{\link{get.real}},\code{\link{fill.covariates}},\code{\link{find.covariates}},\code{\link{inverse.link}},\code{\link{deriv_inverse.link}}
#' @keywords utility
#' @examples
#' 
#' # see examples in fill.covariates
#' 
compute.real <-
function(model,beta=NULL,design=NULL,data=NULL,se=TRUE,vcv=FALSE)
{
# ------------------------------------------------------------------------------------------------
#
#   compute.real  -   computes real estimates and var-cov from design matrix (design) and coefficients
#                     (beta) and inverse link transformation.
#       
# Arguments:  
#
#   model          - MARK model object
#   beta           - estimates of beta parameters for computation of real parameters
#   design         - design matrix from a MARK model
#   data           - dataframe with covariate values that are averaged for estimates
#   se             - if TRUE returns std errors of real estimates
#   vcv            - logical; fct computes and returns v-c matrix of real estimates if TRUE
#
# Value (list):
#
#   real        - data frame containing estimates and se
#   vcv.real    - variance-covariance matrix of real estimates 
#
# ------------------------------------------------------------------------------------------------
model=load.model(model)
if(is.null(beta))beta=model$results$beta$estimate
#
# Assign value of chat - overdispersion constant
#
if(is.null(model$chat))
  chat=1
else
  chat=model$chat
#
# If no design matrix was specified, use the model matrix and fill in any covariate values with default mean value
#
if(is.null(design))
{
   if(is.null(data))
       stop("\n data argument must be specified if design matrix is not specified\n")
   if(!is.data.frame(data))
       stop("\n data argument must be a dataframe. Use data and not processed data list.\n")
   design=fill.covariates(model,find.covariates(model,data))
}
#
# The following shouldn't happen unless the model and design matrices are mixed between models
#
if(dim(design)[2]!=length(beta))
   stop("Mismatch between number of design columns and length of beta")
#
# Change design matrix values to numeric from character and trap to make sure
# the value for design doesn't contain any unfilled covariate entries.
#
if(any(is.na(suppressWarnings(as.numeric(design)))))
   stop("\nInput design matrix must only have numeric values.  Use find.covariates and fill.covariates to fill in covariate values\n")
design=matrix(as.numeric(design),nrow=dim(design)[1])
#
# Set indices for real parameters that have been fixed and at any boundaries
#
if(!is.null(model$fixed))
{
   if(is.null(model$simplify))
      fixedparms=(1:dim(design)[1])%in%model$fixed$index
   else
      fixedparms=(1:dim(design)[1])%in%model$simplify$pim.translation[model$fixed$index]
}
else
   fixedparms=rep(FALSE,dim(design)[1])
boundaryparms=model$results$real$se==0 & !fixedparms
fixedvalues=rep(NA,nrow(design))
fixedvalues[model$simplify$pim.translation[model$fixed$index]]=model$fixed$value
#
#  Compute real parameters; if neither se or vcv then return vector of real parameters
#
   real=convert.link.to.real(design%*%beta,links=model$links,fixed=fixedvalues)
#
#  Set fixed real parameters to their fixed values
#
   real[fixedparms]=fixedvalues[fixedparms]
#  If no se or vcv requested, return result
if(!vcv & !se)return(data.frame(real=real))
#
# Compute vc matrix for real parameters which needs to be modified for
# any mlogit parameters
#
if(length(model$links)==1)
   deriv.real=deriv_inverse.link(real,design,model$links)
else
   deriv.real=t(apply(data.frame(real=real,x=design,links=model$links),1,
         function(x){deriv_inverse.link(as.numeric(x[1]),as.numeric(x[2:(length(x)-1)]),x[length(x)])}))
vcv.real=deriv.real%*%model$results$beta.vcv%*%t(deriv.real)
#
# If vcv=TRUE, compute v-c matrix and std errors of real estimates
# To handle any mlogit parameters compute pseudo-real estimates using log in place of mlogit
#
  ind=grep("mlogit",model$links,ignore.case=TRUE)
  templinks=model$links
  if(length(ind)>0)
  {
    templinks[ind]="log"
    pseudo.real=as.vector(convert.link.to.real(design%*%beta,links=templinks))
#   IF fixed parameters are included in mlogit set, need to recompute the real parameters using 1 in
#   mlogit calculation for every fixed parameter
	if(any(fixedparms[ind]))
	{
		pseudo.real[fixedparms]=0
		pseudo.real[ind][fixedparms[ind]]=exp(pseudo.real[ind][fixedparms[ind]])
		# sums=by(pseudo.real[ind],model$links[ind],sum)
		
		# replacement keeps levels in original order, 
		# otherwise by function creates factor out of character that rearranges
		# the order, only noticeably if >10 levels because factor will order as 
		# 1, 10, 11, 2, 3, 4...
		sums = by(pseudo.real[ind], factor(model$links[ind], levels = unique(model$links[ind])), sum)
		
		sums=sums[match(model$links[ind],names(sums))]
		
		# real[ind]=pseudo.real[ind]/(1+sums[ind])
		# sums is already of correct length, [ind] gives error because it refers
		# to indices outside the size of sums vector
		# I could be wrong, though.
		real[ind]=pseudo.real[ind]/(1+sums)
		
		real[fixedparms]=fixedvalues[fixedparms]
	}
#
#   Compute first derivatives of pseudo-real (for any mlogit parameter)
#   estimates with respect to beta parameters
#
    if(length(templinks)==1)
      deriv.pseudo=deriv_inverse.link(pseudo.real,design,templinks)
   else
      deriv.pseudo=t(apply(data.frame(real=pseudo.real,x=design,links=templinks),1,
            function(x){deriv_inverse.link(as.numeric(x[1]),as.numeric(x[2:(length(x)-1)]),x[length(x)])}))
    deriv.pseudo[fixedparms,]=0
    vcv.pseudo=chat*deriv.pseudo%*%model$results$beta.vcv%*%t(deriv.pseudo)
#
#    Apply chain rule to get variance of real parameters which has mlogits
#    expressed as zi/(1+z1+...zk) where k is number of mlogit components-1 and
#    non-mlogits are expressed as zi.
#    bottom is either 1 for non-mlogits and the sum for mlogits
#    pbottom is partial with respect to zi
#
    if(length(model$links)==1)
        links=rep(model$links,length(pseudo.real))
    else
        links=model$links
    mlogits=outer(links,links,function(x,y)as.numeric(x==y))*as.numeric(substr(links,1,6)=="mlogit"|substr(links,1,6)=="MLogit")
    pbottom=matrix(0,nrow=dim(vcv.pseudo)[1],ncol=dim(vcv.pseudo)[1]) + mlogits
    bottom=diag(nrow=dim(vcv.pseudo)[1])*(1-as.numeric(substr(links,1,6)=="mlogit"|substr(links,1,6)=="MLogit"))+
       mlogits + pbottom*apply(pbottom*pseudo.real,2,sum)
    deriv.pseudo=(diag(nrow=dim(vcv.pseudo)[1])*bottom-pseudo.real*pbottom)/bottom^2
    deriv.pseudo[is.nan(deriv.pseudo)]=0
    vcv.real=deriv.pseudo%*%vcv.pseudo%*%t(deriv.pseudo)
  }
  else
     vcv.real=chat*vcv.real
#
# Compute conf interval taking into account use of logit transform for mlogit
# and any 0-1 link (loglog,cloglog,sin,logit)
#
link.se=suppressWarnings(sqrt(chat*diag(design%*%model$results$beta.vcv%*%t(design))))
link.se[is.na(link.se)]=0
if(length(model$links)==1)
  links=rep(model$links,length(real))
else
  links=model$links
ind=unique(c(grep("mlogit",model$links,ignore.case=TRUE),which(links%in%c("sin","Sin","LogLog","loglog","CLogLog","cloglog"))))
linkse=suppressWarnings(sqrt(diag(vcv.real)[ind])/(real[ind]*(1-real[ind])))
linkse[is.na(linkse)]=0
linkse[is.infinite(linkse)]=0
link.se[ind]=linkse
link.values=design%*%beta
link.values[ind]=suppressWarnings(log(real[ind]/(1-real[ind])))
link.values[ind][abs(real[ind]-1)<1e-7]=100
link.values[ind][abs(real[ind]-0)<1e-7]=-100
links[ind]="logit"
real.lcl=convert.link.to.real(link.values-1.96*link.se,links=links)
real.ucl=convert.link.to.real(link.values+1.96*link.se,links=links)
#
# Set v-c values of fixed parameters to 0
#
vcv.real[fixedparms,]=0
vcv.real[,fixedparms]=0
diag(vcv.real)[diag(vcv.real)<0]=0
se.real=sqrt(diag(vcv.real))
#se.real[is.na(se.real)]=0
fixed=rep("",dim(design)[1])
fixed[fixedparms]="Fixed"
fixed[boundaryparms]="Boundary"
if(vcv)
   return(list(real=real,se.real=se.real,lcl=real.lcl,ucl=real.ucl,fixed=fixed,vcv.real=vcv.real))
else
   return(data.frame(estimate=real,se=se.real,lcl=real.lcl,ucl=real.ucl,fixed=fixed))
}

