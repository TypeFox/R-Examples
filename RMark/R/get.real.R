#' Extract or compute sets of real parameters
#' 
#' Extracts or computes real parameters for a particular type of parameter
#' (parameter) and returns in either a table (dataframe) format or in PIM
#' format.
#' 
#' This function is called by \code{\link{summary.mark}} to extract(from
#' \code{model$results$real}) sets of parameters for display.  But, it can also
#' be useful to compute particular sets of real parameters for output,
#' manipulation or plotting etc.  It is closely related to
#' \code{\link{compute.real}} and it uses that function when it computes
#' (rather than extracts) real parameters. It provides an easier way to
#' extract/compute real estimates of a particular type (\code{parameter}).
#' 
#' The real parameter estimates are computed when either 1) \code{model$chat} >
#' 1, 2) \code{design}, \code{data}, or \code{beta} are specified with
#' non-NULL values for those arguments, or 3) vcv=TRUE.  If none of the above
#' hold, then the estimates are extracted.
#' 
#' If \code{se=FALSE} and estimates are shown in triangular or square PIM
#' structure depending on the model and parameter.  For triangular, the lower
#' half of the matrix is shown as NA (not applicable in this case). If
#' \code{se=TRUE}, the estimate, standard error and confidence interval for the
#' real parameters with the accompanying design data are combined into a
#' dataframe.
#' 
#' If the model contains individual covariates, there are 3 options for
#' specifying the covariate values that are used for the real estimates.  If
#' neither \code{design} nor \code{data} are specified, then the real estimates
#' are computed with the average covariate values used in the MARK output.
#' This is what is done in the call from \code{\link{summary.mark}}.
#' Alternatively, the argument \code{design} can be given a numeric design
#' matrix for the model with the covariates given specific values.  This can be
#' done with \code{\link{find.covariates}} and \code{\link{fill.covariates}}.
#' Finally, a quicker approach is to specify a dataframe for \code{data} which
#' is averaged for the numeric covariate values and these are automatically
#' filled into the design matrix of the model with calls to
#' \code{\link{find.covariates}} and \code{\link{fill.covariates}} from
#' \code{\link{compute.real}} which is used for computation.  The second and
#' third options are essentially the same but creating the completed design
#' matrix will be quicker if multiple calls are made with the same completed
#' design matrix for different parameters.  The dataframe for \code{data} can
#' contain a single entry to specify certain values for computation.
#' 
#' @param model MARK model object
#' @param parameter type of parameter in model (character) (e.g.,"Phi")
#' @param beta values of beta parameters for computation of real parameters
#' @param se if TRUE uses table format and extracts se and confidence
#' intervals, if FALSE uses PIM format with estimates only
#' @param design a numeric design matrix with any covariate values filled in
#' with numerical values
#' @param data covariate data to be averaged for estimates if design=NULL
#' @param vcv if TRUE computes and returns the v-c matrix of the subset of the
#' real parameters
#' @param show.fixed if TRUE fixed values are returned rather than NA in place
#' of fixed values
#' @param expand if TRUE, returns vcv matrix for unique parameters and only 
#' simplified unique parameters if FALSE
#' @param pim if TRUE and se=FALSE, returns a list of estimates in PIM format
#' @return estimates: if \code{se=FALSE and Beta=NULL}, a matrix of estimates
#' or list of matrices for more than one group, and if \code{se=TRUE or beta=is
#' not NULL and vcv=FALSE} a dataframe of estimates with attached design data.
#' If \code{vcv=TRUE}, a list is returned with elements \code{vcv.real} and the
#' dataframe \code{estimates} as returned with \code{se=TRUE}.
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{summary.mark}},\code{\link{compute.real}}
#' @keywords utility
#' @examples
#' 
#' data(example.data)
#' pregion=list(formula=~region)
#' PhiAge=list(formula=~Age)
#' mod=mark(example.data,model.parameters=list(p=pregion,Phi=PhiAge),
#'  groups=c("sex","age","region"),age.var=2,initial.ages=c(0,1,2),threads=1)
#' # extract list of Phi parameter estimates for all groups in PIM format 
#' Phi.estimates=get.real(mod,"Phi")  
#' # print out parameter estimates in triangular PIM format
#' for(i in 1:length(Phi.estimates))  
#' {
#'   cat(names(Phi.estimates)[i],"\n")
#'   print(Phi.estimates[[i]]$pim,na.print="")
#' }
#' require(plotrix)
#' #extract parameter estimates of capture probability p with se and conf intervals
#' p.table=get.real(mod,"p",se=TRUE) 
#' print(p.table[p.table$region==1,])  # print values from region 1
#' estimates=by(p.table$estimate,p.table$region,mean)
#' lcl=by(p.table$lcl,p.table$region,mean)
#' ucl=by(p.table$ucl,p.table$region,mean)
#' plotCI(c(1:4),estimates,ucl-estimates,estimates-lcl,xlab="Region",
#'          ylab="Capture probability",
#' 		ylim=c(.5,1),main="Capture probability estimates by region")
#' 
#' 
get.real <-
function(model,parameter,beta=NULL,se=FALSE,design=NULL,data=NULL,vcv=FALSE,show.fixed=TRUE,expand=FALSE,pim=TRUE)
{
#
# Functions used: valid.parameters, Put.in.PIM.format (defined internally below), compute.real,
#                 setup.parameters
# -----------------------------------------------------------------------------------------------
#
# Internal function to create PIM format
#
  Put.in.PIM.format=function(real,pim,design.data)
  {
    wtable=matrix(NA,nrow=nrow(pim),ncol=ncol(pim))
    for(i in 1:nrow(pim))
    {
	  start=i%%ncol(pim)
	  if(start==0)
		  cols=ncol(pim)
	  else
		  cols=start:ncol(pim)
#	  cols=(i:ncol(pim))
#	  if(i>ncol(pim))cols=((i-ncol(pim)):ncol(pim))
      wtable[i,cols]=real[pim[i,cols]]
    }
    colnames(wtable)=design.data$time[1:ncol(pim)]
    if(!is.null(design.data$cohort))
		if(nrow(pim)==ncol(pim))
           rownames(wtable)=design.data$cohort[(diag(pim[1:nrow(pim),,drop=FALSE])-pim[1,1]+1)]
        else
			rownames(wtable)=rep(design.data$cohort[(diag(pim[1:nrow(pim),,drop=FALSE])-pim[1,1]+1)],nrow(pim)/ncol(pim))	
    wtable
  }
#
# First check to make sure model has been run
#
  model=load.model(model)
  if(is.null(model$results)) 
  {
     message("Model output is not available\n")
     invisible()
  }
#
# Next make sure that requested parameter is appropriate
#
  if(!valid.parameters(model$model,parameter))stop()
#
#  Check to see if there are any covariates used in the model.   If there
#  are covariates used and no design is specified but data is given, then create
#  the design matrix from data
#
  if("covariates" %in% names(model) && !is.null(model$covariates))
    if(is.null(design)& !is.null(data))
       design=fill.covariates(model,find.covariates(model,data=data))
#
#If vcv=TRUE set se=TRUE
#
  if(vcv)se=TRUE
#
#If se=TRUE, warn user if beta.vcv has negative variances
#
  if(se)
    if(any(diag(model$results$beta.vcv)<0))
       warning("\nImproper V-C matrix for beta estimates. Some variances non-positive.\n")
#
#  If beta vector specified, set se and vcv=F
#
   if(!is.null(beta))
   {
	   if(se | vcv )message("Note: with beta specified, having se or vcv set to TRUE is not allowed; setting to FALSE")
	   se=FALSE
	   vcv=FALSE
   }
#
# If design is specified, then compute the real estimates with this design matrix
#
  if(!is.null(design))
  {
     if(!se)
     {
       real=compute.real(model,design=design,beta=beta,se=FALSE)$real
       if(!is.null(model$simplify))
           real=real[model$simplify$pim.translation]
     }
     else
     {
       if(!vcv)
          real=compute.real(model,design=design,beta=beta,vcv=FALSE)
       else
       {
          real.list=compute.real(model,design=design,beta=beta,vcv=TRUE)
          real=data.frame(estimate=real.list$real,se=real.list$se.real,
                          lcl=real.list$lcl,ucl=real.list$ucl)
       }
       if(!is.null(model$simplify))
       {
          real=real[model$simplify$pim.translation,]
          rownames(real)=model$simplify$real.labels
       }
     }
  }
#
# If design and data are not specified, extract estimates (se=FALSE) or
# estimate table (includes se and ci); if model$chat is defined and not 1
# then recompute the confidence intervals for the real parameters with the
# se for beta inflated by sqrt(model$chat). Also if vcv=TRUE recompute to get
# vcv for reals.
#
  else
  {
     if(!se & is.null(beta))
       if(!is.null(model$simplify))
           real=model$results$real$estimate[model$simplify$pim.translation]
       else
           real = model$results$real$estimate
     else
     {
        if(vcv || !is.null(beta) || (!is.null(model$chat) && model$chat !=1))
        {
           if(is.null(model$results$covariate.values$Value))
              data=NULL
          else
              data = as.data.frame(t(model$results$covariate.values$Value))
           if (!is.null(data) && dim(data)[1] != 0)
           {
              names(data)=model$covariates
              if(!vcv)
                 real=compute.real(model,data=data,beta=beta,vcv=FALSE)
              else
              {
                  real.list=compute.real(model,data=data,beta=beta,vcv=TRUE,se=se)
                  real=data.frame(estimate=real.list$real,se=real.list$se.real,
                          lcl=real.list$lcl,ucl=real.list$ucl)
              }
           }
           else
           {
              if(!vcv)
                 real=compute.real(model,design=model$design.matrix,beta=beta,vcv=FALSE,se=se)
              else
              {
                   real.list=compute.real(model,design=model$design.matrix,beta=beta,vcv=TRUE,se=se)
                   real=data.frame(estimate=real.list$real,se=real.list$se.real,
                          lcl=real.list$lcl,ucl=real.list$ucl)
              }
           }
        }
        else
           real=model$results$real
        if(!is.null(model$simplify))
        {
		   if(!se & is.data.frame(real))real=real[,1]
		   if(is.vector(real))
		   {
			   real=real[model$simplify$pim.translation]
			   names(real)=model$simplify$real.labels
		   } else
		   {
			   real=real[model$simplify$pim.translation,,drop=FALSE]
			   rownames(real)=model$simplify$real.labels
		   }
        }
     }
  }
#
#  Set non-estimated reals to NA or fixed value depending on show.fixed
#
   if(!is.null(model$fixed))
	 if(is.vector(real))
     if(show.fixed)
		real[model$fixed$index]=model$fixed$value
     else
        real[model$fixed$index]=NA

  else
     if(is.data.frame(real))
     {
       if(show.fixed)
       {
		   real$estimate[model$fixed$index]=model$fixed$value
		   real$se[model$fixed$index]=0
       }
       else
       {
		   real$estimate[model$fixed$index]=NA
		   real$se[model$fixed$index]=NA
       }
	   real$lcl[model$fixed$index]=real$estimate[model$fixed$index]
	   real$ucl[model$fixed$index]=real$estimate[model$fixed$index]
     }
#
# Lookup type of PIM structure
#
  parameters=setup.parameters(model$model)
  parameter.names=names(parameters)
  type=parameters[[match(parameter,parameter.names)]]$type
  ng=length(model$pims[[parameter]])
  if( type%in%c("Triang","STriang") | !is.null(model$mixtures)| !is.null(model$nocc.secondary))
  {
	  if(!se & pim &is.null(beta))
		  estimates=list()
	  else
		  estimates=data.frame()
  } else
      estimates=NULL
  parameter.labels=names(model$pims[[parameter]][[1]])
  parameter.labels=parameter.labels[!parameter.labels%in%"pim"]
  if(length(model$strata.labels)==1) parameter.labels=parameter.labels[!parameter.labels%in%"stratum"]
  if(model$number.of.groups==1) parameter.labels=parameter.labels[!parameter.labels%in%"group"]
  output.labels=vector(length=ng)
  for(j in 1:ng)
  {
    if(length(parameter.labels)!=0)
    {
        output.labels[j]=""
        for(k in 1:length(parameter.labels))
        {
           if(parameter.labels[k]=="group")
              output.labels[j]=paste(output.labels[j],"Group:",model$group.labels[model$pims[[parameter]][[j]]$group],sep="")
           else
              if(parameter.labels[k]=="stratum")
                 output.labels[j]=paste(output.labels[j]," Stratum:",model$pims[[parameter]][[j]]$stratum,sep="")
#			     output.labels[j]=paste(output.labels[j]," Stratum:",model$strata.labels[model$pims[[parameter]][[j]]$stratum],sep="")
			 else
                 if(parameter.labels[k]=="tostratum")
					 output.labels[j]=paste(output.labels[j]," To:",model$pims[[parameter]][[j]]$tostratum,sep="")
#				 output.labels[j]=paste(output.labels[j]," To:",model$strata.labels[model$pims[[parameter]][[j]]$tostratum],sep="")
              else
                 if(parameter.labels[k]=="session")
                    output.labels[j]=paste(output.labels[j]," Session:",model$pims[[parameter]][[j]]$session.label,sep="")
				
       }
    }
    else
        output.labels[j]=""
    if(!se & is.null(beta) & pim)
    {
        if(type%in%c("Triang","STriang"))
        
        {
           estimates[[j]]=list(pim=Put.in.PIM.format(real,model$pims[[parameter]][[j]]$pim,model$design.data[[parameter]]))
           if(length(parameter.labels)!=0)
              estimates[[j]]=c(estimates[[j]],model$pims[[parameter]][[j]][parameter.labels])
        }
        else
           if(!is.null(model$mixtures) | !is.null(model$nocc.secondary))
           {
              estimates[[j]]=matrix(real[model$pims[[parameter]][[j]]$pim],nrow=dim(model$pims[[parameter]][[j]]$pim)[1])
              if(!is.matrix(estimates[[j]]))estimates[[j]]=matrix(estimates[[j]])
              colnames(estimates[[j]])=model$design.data[[parameter]]$time[model$pims[[parameter]][[j]]$pim[1,]-min(model$pims[[parameter]][[1]]$pim[1,])+1]
              if(is.null(colnames(estimates[[j]])))colnames(estimates[[j]])=rep("",ncol(estimates[[j]]))
              if(!is.null(model$mixtures)&& model$mixtures>1 && !is.null(model$parameters[[parameter]]$rows))
			      rownames(estimates[[j]])=paste("mixture:",1:(model$mixtures+model$parameters[[parameter]]$rows),sep="")
              else
                  rownames(estimates[[j]])=rep("",nrow(estimates[[j]]))
           }
           else
              estimates=rbind(estimates,real[model$pims[[parameter]][[j]]$pim])
    }
    else
    {
    
        if(type%in%c("Triang","STriang")&&model$parameters[[parameter]]$pim.type!="all")
		   if(model$parameters[[parameter]]$pim.type%in%c("time","age"))
              indices=model$pims[[parameter]][[j]]$pim[1,]
           else
              indices=model$pims[[parameter]][[j]]$pim[1,1]
        else
           indices=as.vector(t(model$pims[[parameter]][[j]]$pim))
        indices=indices[indices>0]
        if(!is.null(model$simplify))
           par.index=model$simplify$pim.translation[indices]
        else
           par.index=indices
        real.names=names(real)
		if(is.vector(real))
			estimate.df=cbind(all.index=indices,par.index,estimate=real[indices])
		else
           estimate.df=cbind(all.index=indices,par.index,real[indices,])
        names(estimate.df)=c("all.diff.index","par.index",real.names)
        estimates=rbind(estimates,estimate.df)
    }
  }
  if(!se& is.null(beta) & pim)
  {
     if(type %in%c("Triang","STriang") | is.list(estimates))
        names(estimates)=output.labels
     else
     {
         if(length(model$begin.time)==1)
            colnames(estimates)=model$design.data[[parameter]]$time[model$pims[[parameter]][[1]]$pim-min(model$pims[[parameter]][[1]]$pim)+1]
         else
            colnames(estimates)=1:dim(estimates)[2]
         rownames(estimates)=output.labels
     }
  } else
  {
	  dd=model$design.data[[parameter]]
	  estimates=cbind(estimates,dd[,!names(dd)%in%c("par.index","model.index")])
	  
  }
#
# Return extracted estimates in chosen format
#
if(!vcv)
   return(estimates)
else
{
   if(!is.null(estimates$par.index))
   {
      if(!expand)
		  pindex=sort(unique(estimates$par.index))
	  else
		  pindex=estimates$par.index 
      real.list$vcv.real=real.list$vcv.real[pindex,pindex,drop=FALSE]
      row.names(real.list$vcv.real)=pindex
      colnames(real.list$vcv.real)=pindex
   }
   return(list(estimates=estimates,vcv.real=real.list$vcv.real))
}
}


