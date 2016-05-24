#' Compute sets of link values for real parameters
#' 
#' Computes link values for real parameters for a particular type of parameter
#' (parameter) and returns in a table (dataframe) format.
#' 
#' This function is very similar to \code{\link{get.real}} except that it
#' provides estimates of link values before they are transformed to real
#' estimates using the inverse-link.  Also, the value is always a dataframe for
#' the estimates and design data and optionally a variance-covariance matrix.
#' See \code{\link{get.real}} for further details about the arguments.
#' 
#' @param model MARK model object
#' @param parameter type of parameter in model (character) (e.g.,"Phi")
#' @param beta values of beta parameters for computation of link values
#' @param design a numeric design matrix with any covariate values filled in
#' with numerical values
#' @param data covariate data to be averaged for estimates if design=NULL
#' @param vcv if TRUE computes and returns the v-c matrix of the subset of the
#' link values
#' @return estimates: If \code{vcv=TRUE}, a list is returned with elements
#' \code{vcv.link} and the dataframe \code{estimates}. If \code{vcv=FALSE},
#' only the estimates dataframe is returned which has the same structure as in
#' \code{\link{get.real}}.
#' @author Jeff Laake
#' @seealso \code{\link{compute.link}},\code{\link{get.real}}
#' @keywords utility
get.link <-
function(model,parameter,beta=NULL,design=NULL,data=NULL,vcv=FALSE)
{
# ----------------------------------------------------------------------------------------
#
# get.link - computes link values for real parameters for a particular type of parameter (parameter) and
#            returns in a table format
#
# Arguments:
#
#  model        - MARK model object
#  parameter    - type of parameter in model (character)
#  beta         - estimates of beta parameters for computation
#  design       - a numeric design matrix with any covariate values filled in with numerical values
#  data         - covariate data to be averaged for estimates if design=NULL
#  vcv          - if TRUE computes and returns the v-c matrix of the subset of the link values
#
# Value:
#  estimates    - a list of tables of estimates (se=TRUE) or a matrices of estimates (se=FALSE)
#                 if more than one group. If only one group the result is a single data.frame(table) or a
#                 matrix.
#
# -----------------------------------------------------------------------------------------------
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
# Compute real parameters that were estimated rowsum (of design matrix) value >0
#
    rowsums=apply(model$design.matrix,1,
      function(x){if(all(x=="0"))return(0) else return(1)})
    rowsums=as.numeric(rowsums>0)
    if(!is.null(model$simplify))rowsums= rowsums[model$simplify$pim.translation]
#
#  Check to see if there are any covariates used in the model.   If there
#  are covariates used and no design is specified but data is given, then create
#  the design matrix from data
#
  if("covariates" %in% names(model) && !is.null(model$covariates))
    if(is.null(design)& !is.null(data))
       design=fill.covariates(model,find.covariates(model,data=data))
#
# If design is specified, then compute the link values with this design matrix
#
  if(!is.null(design))
  {
      link.list=compute.link(model,design=design,beta=beta,vcv=TRUE)
      link=link.list$estimates
      if(!is.null(model$simplify))
      {
          link=link[model$simplify$pim.translation,]
          rownames(link)=model$simplify$real.labels
      }
  }
#
# If design and data are not specified, extract data from the values used
#  in the output or use the mode dseign matrix
#
  else
  {
      data=as.data.frame(model$results$covariate.values$Value)
      if(dim(data)[1]!=0)
      {
         names(data)=model$covariates
         link.list=compute.link(model,data=data,beta=beta,vcv=TRUE)
         link=link.list$estimates
      }
      else
      {
         link.list=compute.link(model,design=model$design.matrix,beta=beta,vcv=TRUE)
         link=link.list$estimates
      }
      if(!is.null(model$simplify))
      {
         link=link[model$simplify$pim.translation,]
         rownames(link)=model$simplify$real.labels
      }
  }
#
#  Set non-estimated links to NA
#
  if(any(rowsums==0))
  {
      link$estimate[rowsums==0]=NA
      link$se[rowsums==0]=NA
      link$lcl[rowsums==0]=link$estimate[rowsums==0]
      link$ucl[rowsums==0]=link$estimate[rowsums==0]
  }
  parameters=setup.parameters(model$model)
  parameter.names=names(parameters)
  type=parameters[[match(parameter,parameter.names)]]$type
  ng=length(model$pims[[parameter]])
  if( type %in%c("Triang","STriang") | !is.null(model$mixtures)| !is.null(model$nocc.secondary))
      estimates=list()
  else
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
                 output.labels[j]=paste(output.labels[j]," Stratum:",model$strata.labels[model$pims[[parameter]][[j]]$stratum],sep="")
              else
                 if(parameter.labels[k]=="tostratum")
                    output.labels[j]=paste(output.labels[j]," To:",model$strata.labels[model$pims[[parameter]][[j]]$tostratum],sep="")
              else
                 if(parameter.labels[k]=="session")
                    output.labels[j]=paste(output.labels[j]," Session:",model$pims[[parameter]][[j]]$session.label,sep="")

       }
    }
    else
        output.labels[j]=""
    if(type%in%c("Triang","STriang") && model$parameters[[parameter]]$pim.type!="all")
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
    link.names=names(link)
    estimate.df=cbind(par.index,link[indices,])
    names(estimate.df)=c("par.index",link.names)
    estimates=rbind(estimates,estimate.df)
  }
  estimates=cbind(estimates,model$design.data[[parameter]])
#
# Return extracted estimates in chosen format
#
if(!vcv)
   return(estimates)
else
{
   if(!is.null(estimates$par.index))
   {
      pindex=sort(unique(estimates$par.index))
      link.list$vcv=link.list$vcv[pindex,pindex]
      row.names(link.list$vcv)=pindex
      colnames(link.list$vcv)=pindex
   }
   return(list(estimates=estimates,vcv.link=link.list$vcv))
}
}
