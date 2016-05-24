#' Extract results from MARK output file (internal use)
#' 
#' Extracts the lnl, AICc, npar, beta and real estimates and returns a list of
#' these results for inclusion in the \code{mark} object. The elements
#' \code{beta} and \code{real} are dataframes with fields estimate,se,lcl,ucl.
#' This function was written for internal use and is called by
#' \code{\link{run.mark.model}}. It is documented here for more advanced users
#' that might want to modify the code or adapt for their own use.
#' 
#' 
#' @param out output from MARK analysis (\code{model$output})
#' @param model mark model object
#' @param adjust if TRUE, adjusts number of parameters (npar) to number of
#' columns in design matrix, modifies AIC and records both
#' @param realvcv if TRUE the vcv matrix of the real parameters is extracted
#' and stored in the model results
#' @param vcvfile name of vcv file output
#' @return result: list of extracted output elements
#' \item{lnl}{-2xLog-likelihood} \item{deviance}{Difference between saturated
#' model and lnl} \item{npar}{Number of model parameters}
#' \item{AICc}{Small-sample corrected AIC value using npar and n}
#' \item{npar.unadjusted}{Number of model parameters as reported by MARK if
#' npar was adjusted} \item{AICc.unadjusted}{Small-sample corrected AIC value
#' using npar.unadjusted and n} \item{n}{Effective sample size reported by
#' MARK; used in AICc calculation} \item{beta}{Dataframe of beta parameters
#' with fields: estimate, se, lcl, ucl} \item{real}{Dataframe of real
#' parameters with fields: estimate, se, lcl, ucl}
#' \item{derived.vcv}{variance-covariance matrix for derived parameters if any}
#' \item{covariate.values}{dataframe with fields Variable and Value which are
#' the covariate names and value used for real parameter estimates in the MARK
#' output} \item{singular}{indices of beta parameters that are non-estimable or
#' at a boundary} \item{real.vcv}{variance-covariance matrix for real
#' parameters (simplified) if realvcv=TRUE}
#' @author Jeff Laake
#' @seealso \code{\link{run.mark.model}}
#' @keywords utility
extract.mark.output <-
function(out,model,adjust,realvcv=FALSE,vcvfile)
{
# ----------------------------------------------------------------------------------------
#
#  extract.mark.output     extracts the lnl, AICc, npar, beta and real estimates 
#                          and returns a list of these
#                          beta and real are dataframes with names estimate,se,lcl,ucl
#  Value:
#     result         - list of extracted output elements like npar, n(ESS), deviance etc
#
#
#  Functions used: setup.model, read.mark.binary 
#
# ----------------------------------------------------------------------------------------
  os=R.Version()$os
#
#  Extract basic stats (npar, n, deviance AICc etc)
#
  locate=function(x)
  {
	  loc=regexpr("} = ",out[x])
	  if(length(loc)==0)return(NULL)
	  if(loc[1]==-1)loc=regexpr(" = ",out[x])-1
	  return(loc)
  }
  design.matrix=model$simplify$design.matrix
  links=model$links
  model_def=setup.model(model$model,model$nocc)
  if(model_def$derived)
	  derived_labels=model_def$derived_labels[[1]]
  else
	  derived_labels=NULL
  outfile=tempfile("markxxx",tmpdir=getwd(),fileext=".tmp")
  nreal=dim(design.matrix)[1]
  nbeta=dim(design.matrix)[2]
  x=grep("Effective sample size ",out,ignore.case=TRUE)
  if(length(x)==0)
  {
	  message("MARK did not run properly.  If error message was not shown, re-run MARK with invisible=FALSE")
	  return(NULL)
  }
  n=type.convert(substr(out[x],regexpr("=",out[x])+1,nchar(out[x])))
  x=grep("-2logL {",out,fixed=TRUE)
  if(length(x)==0)
  {
	  message("MARK did not run properly.  If error message was not shown, re-run MARK with invisible=FALSE")
	  return(NULL)
  }
  #  x=grep("-2logL {",out,ignore.case=TRUE, extended=FALSE)
  lnl=type.convert(substr(out[x],locate(x)+4,nchar(out[x])))
  x=grep("Number of Estimated",out,ignore.case=TRUE)
  if(length(x)==0)
  {
	  message("MARK did not run properly.  If error message was not shown, re-run MARK with invisible=FALSE")
	  return(NULL)
  }
  npar=type.convert(substr(out[x],locate(x)+4,nchar(out[x])))
  x=grep("DEVIANCE ",out,ignore.case=TRUE)
  if(length(x)==0)
  {
	  message("MARK did not run properly.  If error message was not shown, re-run MARK with invisible=FALSE")
	  return(NULL)
  }
  deviance=type.convert(substr(out[x],locate(x)+4,nchar(out[x])))[1]
  x = grep("DEVIANCE Degrees of Freedom ", out, ignore.case = TRUE)
  if(length(x)==0)
  {
	  message("MARK did not run properly.  If error message was not shown, re-run MARK with invisible=FALSE")
	  return(NULL)
  }
  deviance.df = type.convert(substr(out[x], locate(x)+4, nchar(out[x])))[1]
  x=grep("AICc",out,ignore.case=TRUE)
  if(length(x)==0)
  {
	  message("MARK did not run properly.  If error message was not shown, re-run MARK with invisible=FALSE")
	  return(NULL)
  }
  AICc=type.convert(substr(out[x],locate(x)+4,nchar(out[x])))
  if(length(links)==1)
     x1=grep(paste(links,"link"),out,ignore.case=TRUE)
  else 
     x1=grep("parm-specific link",out,ignore.case=TRUE)
 if(length(x1)==0)
 {
	 message("MARK did not run properly.  If error message was not shown, re-run MARK with invisible=FALSE")
	 return(NULL)
 }
 x2=grep("Real Function Parameters",out,ignore.case=TRUE)
  if(length(x2)==0)
  {
	  message("MARK did not run properly.  If error message was not shown, re-run MARK with invisible=FALSE")
	  return(NULL)
  }
  x3=length(out)
  if(length(grep("proc stop",out,ignore.case=TRUE))==0)
     message("\nWarning: output from MARK was not complete\n")
  x4=grep("Variable   Value",out,ignore.case=FALSE)+1
  if(length(x4)==0)x4=x2
#
# Extract average covariate values used in real parameter calculation
#
  if(x4>x2)
  {
     ff <- tempfile()
     cat(file=ff, out[(x4+1):(x4+length(model$covariates))],sep="\n")
     covariate.values=read.fwf(file=ff,widths=c(20,15),col.names=c("Variable","Value"))
  }
  else
     covariate.values=NULL
#
# Extract beta parameters ; this could also be done from binary file but from text
# file it is easy to decide which values are fixed.
#
  j=1
  save=NULL
  for(i in x1:(x2-1))
  {
    if(j<=nbeta)
    {
       ind=regexpr(paste(" ",j,":",sep=""),out[i])
       if(ind!=-1 & ind<=20)
       {
         save=c(save,out[i])
         j=j+1
       }
     }
  }
  write(save,outfile)
  x=read.fwf(file=outfile,widths=c(26,16,16,16,16),col.names=c("","estimate","se","lcl","ucl"))
  dimx=dim(x)[2]
  beta=as.data.frame(x[,((dimx-4+1):dimx)])
  names(beta)=c("estimate","se","lcl","ucl")
  row.names(beta)= colnames(model$simplify$design.matrix)
  nbeta=length(beta$estimate[beta$estimate!=0.000000])
#
# Extract parameter numbers that were not "estimated"
#
  singular=NULL
  if(nbeta!=npar)
  {
     x=grep("Attempted ordering of parameters",out,ignore.case=TRUE)
     if(length(x)==0)
       warning("\nNot all parameters were estimated but not able to find non-estimable parameters\n")
     else
     {
        nlines=ceiling(nbeta/25)
        par.indices=NULL
        for (i in (x+1):(x+nlines))
		{
			ii=strsplit(out[i]," ")[[1]]
			par.indices=c(par.indices,as.numeric(ii[ii!=""]))
		}
        singular=par.indices[(npar+1):nbeta]
     }     
  }
  if(nbeta!=npar & adjust)
  {
    message("\nNote: only ",npar," parameters counted of ",nbeta," specified parameters\n") 
    message("AICc and parameter count have been adjusted upward\n")
    AICc.unadjusted=AICc
    npar.unadjusted=npar    
    AICc=lnl+ 2*nbeta +2*nbeta*(nbeta+1)/(n - nbeta -1)
    npar=nbeta
  }
  else
    npar.unadjusted=NULL    
  unlink(outfile)
#
# Extract real parameters from text file; This could be done from binary file but the text file
# also denotes the fixed parameters
#
  j=1
  if(x4>x2)x2=x4+length(model$covariates)+1
  for(i in x2:(x3-1))
  {
    if(j<=nreal)
    {
       ind=regexpr(paste(" ",j,":",sep=""),out[i])
       if(ind==-1)
		   ind= regexpr("\\*\\*\\*\\*:",out[i])
       if(ind!=-1& ind<=20)
       {
          write(out[i],file=outfile,append=TRUE)
          j=j+1
       }
    }
  }

  x=read.fwf(file=outfile,widths=c(27,16,16,16,16,20),col.names=c("","estimate","se","lcl","ucl","fixed"),
                               as.is=TRUE)
  unlink(outfile)
  x$note=""
  x$fixed[is.na(x$fixed)]=  "       "
  x$note[substr(as.character(x$fixed),3,7)!="Fixed"]=x$fixed[substr(as.character(x$fixed),3,7)!="Fixed"]
  x$fixed[substr(as.character(x$fixed),3,7)!="Fixed"]="     "
  x$fixed[substr(as.character(x$fixed),3,7)=="Fixed"]="Fixed"
  x$fixed[is.na(x$fixed)]=  "       "
  real=data.frame(estimate=as.numeric(x$estimate),se=as.numeric(x$se),lcl=as.numeric(x$lcl),ucl=as.numeric(x$ucl),fixed=x$fixed,note=x$note)
  if(is.null(model$simplify))
     row.names(real) = row.names(design.matrix)
  else
     row.names(real)=row.names(model$simplify$design.matrix)
#  if(!is.factor(real$fixed))real$fixed=""
  if(!is.factor(real$note))real$note=""
  if(file.exists(vcvfile))
     if(os=="mingw32")
        param=read.mark.binary(vcvfile,derived_labels)
     else
        param=read.mark.binary.linux(vcvfile,derived_labels)
  else
  {
     param=NULL
     message("\nV-C file is missing. Skipping over it.\n")
  }
  if(realvcv)
    real.vcv=param$real.vcv
  else
    real.vcv=NULL
  if(is.null(npar.unadjusted))
     return(list(lnl=lnl,deviance=deviance,deviance.df=deviance.df,npar=npar,n=n,AICc=AICc,beta=beta,real=real,beta.vcv=param$beta.vcv,derived=param$derived,derived.vcv=param$derived.vcv,
                 covariate.values=covariate.values,singular=singular,real.vcv=real.vcv))

  else
     return(list(lnl=lnl,deviance=deviance,deviance.df=deviance.df,npar=npar,npar.unadjusted=npar.unadjusted,n=n,AICc=AICc,AICc.unadjusted=AICc.unadjusted,
                 beta=beta,real=real,beta.vcv=param$beta.vcv,derived=param$derived,derived.vcv=param$derived.vcv,
                 covariate.values=covariate.values,singular=singular,real.vcv=real.vcv))
}
