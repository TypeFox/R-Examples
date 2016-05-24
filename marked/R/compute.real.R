#' Compute estimates of real parameters
#' 
#' Computes real estimates and their var-cov for a particular parameter. 
#' 
#' This code is complicated because it handles both the MCMC models and the likelihood models. The former is quite 
#' simple than the latter because all of the real computation is done by the model code and this function only computes 
#' summaries.  The likelihood model code is complicated primarily by the mlogit parameters which are computed in 2 stages:
#' 1) log link and 2) summation to normalize. The mlogit is handled differently depending on the model.
#' For MS and JS models, one of the parameters is computed by subtraction (specified as addone==TRUE) whereas the
#' HMM models (addone=FALSE) specify a parameter for each cell and one is fixed by the user to 1. The latter is preferable
#' because it then provides an estimate and a std error for each parameter whereas the subtracted value is not provided for MS and JS.  
#' 
#' This function differs from compute.real in RMark because it only computes the values for a single parameter whereas the function 
#' with the same name in RMark can compute estimates from multiple parameters (eg Phi and p). 
#' 
#' @param model model object
#' @param parameter name of real parameter to be computed (eg "Phi" or "p")
#' @param ddl list of design data 
#' @param dml design matrix list
#' @param unique TRUE if only unique values should be returned
#' @param vcv logical; if TRUE, computes and returns v-c matrix of real estimates
#' @param se logical; if TRUE, computes std errors and conf itervals of real estimates
#' @param chat over-dispersion value
#' @param subset logical expression using fields in real dataframe
#' @param select character vector of field names in real that you want to include 
#' @param showDesign if TRUE, show design matrix instead of data 
#' @param include vector of field names always to be included even when select or unique specified
#' @param uselink default FALSE; if TRUE uses link values in evaluating uniqueness
#' @export
#' @return A data frame (\code{real}) is returned if \code{vcv=FALSE};
#' otherwise, a list is returned also containing vcv.real: \item{real}{ data
#' frame containing estimates, and if vcv=TRUE it also contains
#' standard errors and confidence intervals} \item{vcv.real}{variance-covariance matrix of
#' real estimates}
#' @author Jeff Laake
#' @examples
#' data(dipper)
#' dipper.proc=process.data(dipper,model="cjs",begin.time=1)
#' dipper.ddl=make.design.data(dipper.proc)
#' mod.Phisex.pdot=crm(dipper.proc,dipper.ddl,
#'  model.parameters=list(Phi=list(formula=~sex+time),p=list(formula=~1)),hessian=TRUE)
#' xx=compute.real(mod.Phisex.pdot,"Phi",unique=TRUE,vcv=TRUE)
#' @keywords utility
compute.real <-function(model,parameter,ddl=NULL,dml=NULL,unique=TRUE,vcv=FALSE,se=FALSE,chat=1,subset=NULL,select=NULL,showDesign=FALSE,include=NULL,uselink=FALSE)
{
#  Note that the vector indices has 3 different meanings in the code as the code progresses.
# if ddl not specified return results stored in model
  if(is.null(ddl))return(model$results$reals)
# set mcmc value
  mcmc=ifelse(class(model)[2]=="mcmc",TRUE,FALSE)
# Setup data (ddl for a parameter) and extract dm; currently is either passed as an argument or stored in model
  data=ddl[[parameter]]
  if(!is.null(dml))
	  design=dml[[parameter]]$fe  
  else
	design=model$results$model_data[[paste(parameter,"dm",sep=".")]]
# remove unneeded reals - values that occurred before cohort was released; here indces represents the real
# parameters to be used
  if(!is.null(data$Time)&model$model.parameter[[parameter]]$type=="Triang")
  {
	  indices=which(data$Time>=data$Cohort)
	  data=data[indices,,drop=FALSE]
	  design=design[indices,,drop=FALSE]
  }
# for abundance parameters only use those with non-zero freq
  if(parameter!="N" &!is.null(data$freq))
  {
	  indices=which(data$freq!=0)
	  data=data[indices,,drop=FALSE]
	  design=design[indices,,drop=FALSE]
  }
# design is a sparse matrix and this changes it to a std matrix
  design=as.matrix(design)
# Set up links; for HMM models ulink (utlimate link) is used for mlogit parameters which have a temp "log" link
# ulink is known because there are mlogit variables named in parameters.txt for the parameter. Parameters such as pent in JS and
# Psi in MSCJS, the link is mlogit (rather than log) 
  link=model$model.parameter[[parameter]]$link
  ulink=NULL
  if(link=="mlogit")
  {
	  link=paste("mlogit",apply(data[,model$model.parameter[[parameter]]$mlogit,drop=FALSE],1,paste,collapse=""),sep="")
      addone=TRUE
  } else {
    if(!is.null(model$model.parameter[[parameter]]$mlogit))
	{
		addone=FALSE
		ulink=paste("mlogit",apply(data[,model$model.parameter[[parameter]]$mlogit,drop=FALSE],1,paste,collapse=""),sep="")
	}
  }
  if(length(link)>1)link=link[indices] 
# If there are any fixed parameters (fix in ddl), create boolean for fixed values (fixedparms) and store values (fixedvalues)
  if(!is.null(data$fix))
  {
	  fixedparms=!is.na(data$fix)
	  fixedvalues=data$fix
  }else
  {
	  fixedparms=rep(FALSE,nrow(data))
	  fixedvalues=rep(NA,nrow(data))
  }
# Next select the variables to be extracted for the real parameter calculations; by default this is the variables in the 
# model formula and if any variables for this parameter are needed for mlogit parameter setup (include) these are also selected.
  if(!is.null(ulink))
	  df=cbind(data,fixed=fixedvalues,link=ulink)
  else
	  df=cbind(data,fixed=fixedvalues)
  if(unique)
  {
	  if(!is.null(ulink)&uselink)
  	     varnames=c(select,include,all.vars(model$model.parameters[[parameter]]$formula),"occ","fixed","link")
	 else
 	     varnames=c(select,include,all.vars(model$model.parameters[[parameter]]$formula),"occ","fixed")
  } else
  {
	  if(!is.null(ulink)&uselink)	
  	      varnames=c(include,select,"fixed","link")
	  else
	      varnames=c(include,select,"fixed")
  }
  varnames=unique(varnames)
  if(any(!select%in%names(df))) 
	  warning(paste("These variable names not in data for real estimates: ",paste(varnames[!select%in%names(df)],collapse=""),sep=""))
  varnames=varnames[varnames%in%names(df)]
  df=df[,varnames,drop=FALSE]
# Check to make sure dimensions of beta and design matrix match
  results=model$results
  beta=results$beta[[parameter]]
  if(mcmc)
  {
	  if(ncol(design)!=ncol(model$results$beta.mcmc[[parameter]]))
   	     stop("Mismatch between number of design columns and length of beta")
  } else
  {
	  if(ncol(design)!=length(beta))
         stop("Mismatch between number of design columns and length of beta")
  }
# Compute real parameters using function convert.link.to.real
# added code to handle all parameters being fixed
  if(length(beta)==0)
	  real=fixedvalues
  else
  {
	  if(mcmc)
		  real=convert.link.to.real(t(design%*%t(results$beta.mcmc[[parameter]])),links=link,fixed=fixedvalues)
	  else
		  real=convert.link.to.real(design%*%beta,links=link,fixed=fixedvalues)
     # Set fixed real parameters to their fixed values
     real[fixedparms]=fixedvalues[fixedparms]
  }
# Summarize reals for mcmc models or create dataframe with estimate=real for non-mcmc models
  if(mcmc)
  {
    fitted <- mcmc(real)
    summ <- summary(fitted)
    hpd <- HPDinterval(fitted)
    reals <- data.frame(
		  mode = apply(real, 2, mcmc_mode), 
  	      mean=apply(real, 2, mean), 
          sd = apply(real, 2, sd),
		  CI.lower=hpd[,1], CI.upper=hpd[,2])
  } else 
  {	
	reals=data.frame(estimate=real)
	rownames(reals)=NULL
  }
  # If subset argument is missing (default) skip all subsetting
  if(!is.null(subset))
  {
#   If subset non-blank use it to subset rows
	if(subset!="")
	{
		if(substitute(subset)!="substitute(subset)")
			sub=substitute(subset)
		else
			sub=subset
		indices = eval(sub, data)
    }
	
#   If subset is blank, then if unique is TRUE, use only rows where all values in
#   df are unique.  If not unique, also use all rows.
  }	else
		if(unique)
		{
			if(is.null(select)|| !"occ"%in%select)
   			   indices=which(!duplicated(apply(df[,!names(df)%in%"occ",drop=FALSE],1,paste,collapse="")))
			else
				indices=which(!duplicated(apply(df,1,paste,collapse="")))
		}
		else
			indices = TRUE
# modify indices to make sure and include complete set in mlogit
  if(!is.null(ulink))
  {
	   tlink=unique(ulink[indices])
	   indices=which(ulink%in%tlink)
  } 	   
# subset df, design, link, ulink, reals, real, fixedparms and fixedvalues with the 
# indices from the subsetted rows.
  df=df[indices,,drop=FALSE]	  
  design=design[indices,,drop=FALSE]
  if(length(link)>1)link=link[indices]
  if(ncol(df)==0)df$Intercept=1
  reals=reals[indices,,drop=FALSE]
  real=real[indices]
  fixedparms=fixedparms[indices]
  fixedvalues=fixedvalues[indices]
#  normalize mlogit parameters
  if(!is.null(ulink))
  {
	  ulink=ulink[indices]
	  sums=by(real,ulink,sum)
	  sums=sums[match(ulink,names(sums))]
	  reals=real/sums
  }	  
  df=df[,!names(df)%in%"fixed",drop=FALSE]
  if(showDesign)
	  reals=cbind(design,reals)
  else
	  reals=cbind(df,reals)
#  If vcv or se do computations. for non-mlogit parameters this is accomplished with the
#  calls to compute deriv.real and vcv.real. Below, the confidence interval for each parameter
#  is approximated.
  if(vcv | se)
  {
	 if(is.null(results$beta.vcv))stop("vcv matrix not available in model")
#    Compute vc matrix for real parameters which needs to be modified for
#    any mlogit parameters; for HMM models the link for mlogit parameters is "log"
#    while fitting and for the compuations below until ulink is used.
	 indices=grep(paste(parameter,"\\.",sep=""),colnames(results$beta.vcv))
	 if(substr(link[1],1,6)!="mlogit")
	 {
	    deriv.real=as.matrix(deriv_inverse.link(real,design,link))
	    vcv.real=deriv.real%*%results$beta.vcv[indices,indices,drop=FALSE]%*%t(deriv.real)
	 } 
#   To handle any mlogit parameters compute pseudo-real estimates using log in place of mlogit;
#   This is not needed for HMM mlogit parameters because it is done above
	 if(substr(link[1],1,6)=="mlogit"| !is.null(model$model.parameter[[parameter]]$mlogit))
	 {
	    if(addone)
	    {
		   pseudo.real=as.vector(convert.link.to.real(design%*%beta,links="log"))
#          Compute first derivatives of pseudo-real (for any mlogit parameter)
#          estimates with respect to beta parameters
	       deriv.pseudo=deriv_inverse.link(pseudo.real,design,"log")
	       deriv.pseudo[fixedparms,]=0
	       vcv.pseudo=deriv.pseudo%*%results$beta.vcv[indices,indices,drop=FALSE]%*%t(deriv.pseudo)
		   links=link
	    } else
#         For HMM mlogit parameters, rename deriv.real and vcv.real and compute
#         the real parameters by normalizing so they sum 1.
	    {
		   deriv.pseudo=deriv.real
		   vcv.pseudo=vcv.real
		   pseudo.real=real
		   sums=by(real,ulink,sum)
		   sums=sums[match(ulink,names(sums))]
		   real=real/sums
		   reals$estimate=real
		   links=ulink
		   fixedparms[fixedvalues==1]=FALSE
	    }
#       Apply chain rule to get variance of real parameters which has mlogits
#       expressed as zi/(1+z1+...zk) where k is number of mlogit components-1 
#       pbottom is partial with respect to zi
		pbottom=outer(links,links,function(x,y)as.numeric(x==y))*as.numeric(substr(links,1,6)=="mlogit")
		bottom=as.numeric(addone)*pbottom+pbottom*apply(pbottom*pseudo.real,2,sum)
		deriv.pseudo=(diag(nrow=nrow(vcv.pseudo))*bottom-pseudo.real*pbottom)/bottom^2
		deriv.pseudo[is.nan(deriv.pseudo)]=0
		vcv.real=deriv.pseudo%*%vcv.pseudo%*%t(deriv.pseudo)
	  }
	  vcv.real=chat*vcv.real
#     Compute conf interval taking into account use of logit transform for mlogit
#     and any 0-1 link (loglog,cloglog,sin,logit)
	  link.se=suppressWarnings(sqrt(chat*diag(design%*%results$beta.vcv[indices,indices]%*%t(design))))
	  link.se[is.na(link.se)]=0
	  if(is.null(ulink))
		  if(length(link)==1)
	        links=rep(link,length(real))
	      else
		    links=link
#     Set ind which is index to non logit link (except long link) parameters	 
	  ind=unique(c(grep("mlogit",links,ignore.case=TRUE),which(links%in%c("sin","Sin","LogLog","loglog","CLogLog","cloglog"))))
#     For the set of parmeters indexed by ind transform backwards to logit link from real
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
#     Set v-c values of fixed parameters to 0
	  vcv.real[fixedparms,]=0
	  vcv.real[,fixedparms]=0
	  diag(vcv.real)[diag(vcv.real)<0]=0
	  se.real=sqrt(diag(vcv.real))
	  fixed=rep("",dim(design)[1])
	  fixed[fixedparms]="Fixed"
	  reals=cbind(reals,data.frame(se=se.real,lcl=real.lcl,ucl=real.ucl))   
    }   
#  Setup subset or get unique records 
	#reals=reals[do.call(order, reals),]
	rownames(reals)=NULL
	if(vcv)
		return(list(real=reals,vcv=vcv.real))
    else
        return(reals)
}
