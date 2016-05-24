#' Compute link values from real parameters
#' 
#' Computes link values from reals using 1-1 real to beta(=link)
#' transformation. Also, creates a v-c matrix for the link values if vcv.real
#' is specified.
#' 
#' It has 2 uses both related to model averaged estimates. Firstly, it is used
#' to transform model averaged estimates so the normal confidence interval can
#' be constructed on the link values and then back-transformed to real space.
#' The second function is to enable parametric bootstrapping in which the error
#' distbution is assumed to be multivariate normal for the link values. From a
#' single model, the link values are easily constructed from the betas and
#' design matrix so this function is not needed.  But for model averaging there
#' is no equivalent because the real parameters are averaged over a variety of
#' models with the same real parameter structure but differing design
#' structures.  This function allows for link values and their var-cov matrix
#' to be created from the model averaged real estimates.
#' 
#' @param x vector of real estimates to be converted to link values
#' @param model MARK model object used only to obtain model structure/links
#' etc.  If function is being called for model averaged estimates, then any
#' model in the model list used to construct the estimates is sufficient
#' @param parm.indices index numbers from PIMS for rows in design
#' matrix(non-simplified indices); x[parm.indices] are computed
#' @param vcv.real v-c matrix for the real parameters
#' @param use.mlogits logical; if FALSE then parameters with mlogit links are
#' transformed with logit rather than mlogit for creating confidence intervals
#' for each value
#' @return A list with the estimates (link values) and the links that were
#' used.  If vcv.real = TRUE, then the v-c matrix of the links is also
#' returned.
#' @author Jeff Laake
#' @seealso \code{\link{model.average}}
#' @keywords utility
compute.links.from.reals=function(x,model,parm.indices=NULL,vcv.real=NULL,use.mlogits=TRUE)
{
#
# Computes link values from reals using 1-1 real to beta(=link) transformation. Also,
# creates a v-c matrix for the link values if vcv.real is specified.
#
#  Arguments:
#
#   x     - vector of real estimates
#   model - MARK model object used only to obtain model structure/links etc.  If
#           function is being called for model averaged estimates, then amy model in
#           the model list used to construct the estimates is sufficient.
#
#  parm.indices - non-simplified indices for the real parameters in the model structure
#
#  vcv.real - v-c matrix for the real parameters
#
#  use.mlogits - logical; if FALSE then parameters with mlogit links are transformed
#                with logit rather than mlogit for creating conf intervals for each value
#
#
#  Value:
#
#    A list with the estimates (link values) and  the links that were used.  If
#    vcv.real = TRUE, then the v-c matrix of the links is also returned.
#
# Internal function used to compute links from reals
#
flink=function(x,link)
{
switch(link,
logit=log(x/(1-x)),
log=log(x),
loglog=-log(-log(x)),
cloglog=-log(-log(1-x)),
identity=x,
sin=asin(2*x-1),
Logit=log(x/(1-x)),
Log=log(x),
LogLog=-log(-log(x)),
CLogLog=-log(-log(1-x)),
Identity=x,
Sin=asin(2*x-1)
)
}
#
# Internal function used to compute derivatives of links with respect to reals for
# use in the delta method to construct the v-c for links from the v-c of reals.
#
deriv.link=function(real,link)
{
real=as.vector(real)
switch(link,
logit=1/real+1/(1-real),
log=1/real,
loglog=1/(real*log(real)),
cloglog=-1/(log(1-real)*(1-real)),
identity=1,
sin=1/sqrt(real-real^2),
Logit=1/real+1/(1-real),
Log=1/real,
LogLog=1/(real*log(real)),
CLogLog=-1/(log(1-real)*(1-real)),
Identity=1,
Sin=1/sqrt(real-real^2)
)
}
#
# Set up links; for mlogit parameters the 95% ci is set with logit transform
# and this is determined by argument use.mlogits. It should be FALSE to set ci on
# mlogits but should be TRUE to get actual betas and their v-c say for simulation
#
model=load.model(model)
links=model$links
if(!is.null(model$simplify$links))links=model$simplify$links
#
# To mimic MARK, a logit transform is used for any link that constrains
# estimates to be probabilities. (sin,logit,cloglog,loglog)
#
links[which(model$links%in%c("logit","Logit","sin","Sin","LogLog","loglog","CLogLog","cloglog"))]="logit"
if(!use.mlogits)
{
   ind=grep("mlogit",links,ignore.case=TRUE)
   links[ind]="logit"
}
if(!is.null(parm.indices))
   if(length(links)>1)
      links=links[parm.indices]
#
#  Fixed parameters will have Nan entries in v-c matrix. These are set to 0
#
vcv.real[is.nan(vcv.real)]=0
#
# Compute link values: The first is the case in which all params use same link function
#
   if(length(links)==1)
   {
      link.values=flink(x,links)
#
#     If vcv.real is given then compute the v-c for links.  The derivative matrix is
#     diagonal in this case because each link is only a function of one real.
#
      if(!is.null(vcv.real))
      {
         deriv.link.values=deriv.link(x,links)
         deriv.link.matrix=matrix(0,nrow=length(deriv.link.values),ncol=length(deriv.link.values))
         diag(deriv.link.matrix)=deriv.link.values
         vcv.links=deriv.link.matrix%*%vcv.real%*%t(deriv.link.matrix)
      }
   }
   else
#
# If they are not all the same but none use multinomial logit (mlogit) then use apply to use
# different link for each parameter.  This is similar to above in that the derivative matrix is
# diagonal.
#
   {
      ind=grep("mlogit",links,ignore.case=TRUE)
      if(length(ind)==0)
      {
         link.values=apply(data.frame(x=x,links=links),1,function(x){flink(as.numeric(x[1]),x[2])})
         if(!is.null(vcv.real))
         {
            deriv.link.values=apply(data.frame(x=x,links=links),1,function(x){deriv.link(as.numeric(x[1]),x[2])})
            deriv.link.matrix=matrix(0,nrow=length(deriv.link.values),ncol=length(deriv.link.values))
            diag(deriv.link.matrix)=deriv.link.values
            vcv.links=deriv.link.matrix%*%vcv.real%*%t(deriv.link.matrix)
         }
      }
      else
#
# If some use multinomial logit (mlogit) then the inverses for those parameters use the log link and then
# they are summed and the mlogit is computed as exp()/(1+sum(exp()).  The link is then log(real/(1-sum(reals)).
# For the other parameters using different links they are applied directly. The variable "ind" indexes those that
# contain mlogit.  The derivative matrix is not diagonal for the portions of the matrix involving mlogit parameters
# because each link computation is composed of all the reals asscociated with that mlogit.  The addon calculation
# copes with those non-diagonal terms.
#
      {
        newlinks=links
        newlinks[ind]="log"
        link.values=apply(data.frame(x=x,links=newlinks),1,function(x){flink(as.numeric(x[1]),x[2])})
        sums=by(x,links,sum)
        sums=sums[match(links,names(sums))]
        link.values[ind]=link.values[ind]-log(1-sums[ind])
        if(!is.null(vcv.real))
        {
           deriv.link.values=apply(data.frame(x=x,links=newlinks),1,function(x){deriv.link(as.numeric(x[1]),x[2])})
           deriv.link.matrix=matrix(0,nrow=length(deriv.link.values),ncol=length(deriv.link.values))
           diag(deriv.link.matrix)=deriv.link.values
           addon=rep(0,length(deriv.link.values))
           addon[ind]=1/(1-sums[ind])
           deriv.link.matrix=deriv.link.matrix+addon*outer(links,links,"==")
           vcv.links=deriv.link.matrix%*%vcv.real%*%t(deriv.link.matrix)
      }
      }
   }
#
#  If parm.indices were specified use them as names
#
   if(!is.null(parm.indices))
   {
      names(link.values)=parm.indices
      if(!is.null(vcv.real))
      {
         row.names(vcv.links)=parm.indices
         colnames(vcv.links)=parm.indices
      }
   }
#  return values and optionally vcv for links
   if(is.null(vcv.real))
      return(list(estimates=link.values,links=links))
   else
      return(list(estimates=link.values,links=links,vcv=vcv.links))
}


