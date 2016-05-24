#' Convert link values to real parameters
#' 
#' Computes real parameters from link values
#' 
#' Computation of the real parameter from the link value is relatively
#' straightforward for most links and the function \code{\link{inverse.link}}
#' is used.  The only exception is parameters that use the \code{mlogit} link
#' which requires the transformation across sets of parameters.  This is a
#' convenience function that does the necessary work to convert from link to
#' real for any set of parameters.  The appropriate links are obtained from
#' \code{model$links} unless the argument \code{links} is specified and they
#' will over-ride those in \code{model}.
#' 
#' @param x Link values to be converted to real parameters
#' @param model model object
#' @param links vector of character strings specifying links to use in
#' computation of reals
#' @param fixed vector of fixed values for real parameters that are needed for
#' calculation of reals from mlogits when some are fixed
#' @return vector of real parameter values
#' @author Jeff Laake
#' @seealso \code{\link{inverse.link}},\code{\link{compute.real}}
#' @keywords utility
convert.link.to.real=function(x,model=NULL,links=NULL,fixed=NULL)
{
#
# Arguments
#  x      Link values to be converted to real parameters
#  model  MARK model object
#  links  vector of character strings specifying links to use in computation of reals
#
# Value:
#  vector of real parameter estimates
#
# Computation of the real parameter from the link value is relatively straightforward for most links and
# the function inverse.link is used.  The only exception is parameters that use the
#  mlogitlink which requires the transformation across sets of parameters.  This is a
# convenience function that does the necessary work to convert from link to real for any set of
# parameters.  The appropriate links are obtained from model$links unless the argument
# links is specified and they will over-ride those in model.
#
if(is.null(links))links=model$links
#                                           ((
# Compute real values: The first is the case in which all params use same link function
#
   if(length(links)==1)
      real=inverse.link(x,links)
   else
#
# If they are not all the same but none use multinomial logit (mlogit) then use apply to use
# different link for each parameter.
#
   {
      ind=grep("mlogit",links,ignore.case=TRUE)
      if(length(ind)==0)
         real=apply(data.frame(x=x,links=links),1,function(x){inverse.link(as.numeric(x[1]),x[2])})
      else
#
# If some use multinomial logit (mlogit) then the inverses for those parameters use the log link and then
# they are summed and the mlogit is computed as exp()/(1+sum(exp()).  For the other parameters using different
# links they are applied directly. The variable "ind" indexes those that contain mlogit.
#
      {
        newlinks=links
        newlinks[ind]="log"
        x[!is.na(fixed)]=fixed[!is.na(fixed)]
        real=apply(data.frame(x=x,links=newlinks),1,function(x){inverse.link(as.numeric(x[1]),x[2])})
        sums=by(real,links,sum)
        sums=sums[match(links,names(sums))]
        real[ind]=real[ind]/(1+sums[ind])
      }
   }
#
# Reset values for fixed parameters
#
   if(!is.null(model))
      real[model$results$real$se==0]=model$results$real$estimate[model$results$real$se==0]
   return(real)
}
