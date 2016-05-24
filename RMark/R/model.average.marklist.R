#' Compute model averaged estimates of real parameters
#' 
#' Computes model averaged estimates and standard errors of real parameters for
#' a list of models with a \code{model.table} constructed from
#' \code{\link{collect.models}}. It can also optionally compute the var-cov
#' matrix of the averaged parameters and their confidence intervals by
#' transforming with the link functions, setting normal confidence intervals on
#' the transformed values and then back-transforming for the real estimates.
#' 
#' If there are any models in the \code{model.list} which do not have any
#' output or results they are dropped.  If any have non-positive variances for
#' the betas and \code{drop=TRUE}, then the model is reported and dropped from
#' the model averaging.  The weights are renormalized for the remaining models
#' that are not dropped before they are averaged.
#' 
#' If \code{parameter=NULL}, all real parameters are model averaged but the
#' design data is not copied over because it can vary by the type of parameter.
#' It is only necessary to model average all parameters at once to get
#' covariances of model averaged parameters of differing types.
#' 
#' If \code{data=NULL}, the average covariate values are used for any models
#' using covariates. Note that this will only work with models created after
#' v1.5.0 such that average covariate values are stored in each model object.
#' 
#' @usage \method{model.average}{marklist}(x, parameter, data, vcv, drop=TRUE, indices=NULL, revised=TRUE, mata=FALSE,
#'         normal.lm=FALSE, residual.dfs=0, alpha=0.025,...)
#' @export 
#' @param x a list of mark model results and a model.table constructed by
#' \code{\link{collect.models}}
#' @param parameter name of model parameter (e.g., "Phi" for CJS models); if
#' left NULL all real parameters are averaged
#' @param data dataframe with covariate values that are averaged for estimates
#' @param vcv logical; if TRUE then the var-cov matrix and confidence intervals
#' are computed
#' @param drop if TRUE, models with any non-positive variance for betas are
#' dropped
#' @param indices a vector of parameter indices from the all-different PIM
#' formulation of the parameter estimates that should be presented.  This
#' argument only works if the parameter argument = NULL.  The primary purpose
#' of the argument is to trim the list of parameters in computing a vcv matrix
#' of the real parameters which can get too big to be computed with the
#' available memory
#' @param revised if TRUE, uses revised variance formula (eq 6.12 from Burnham
#' and Anderson) for model averaged estimates and eq 6.11 when FALSE
#' @param mata if TRUE, create model averaged tail area confidence intervals as described by Turek and Fletcher
#' @param alpha The desired lower and upper error rate.  Specifying alpha=0.025
#' corresponds to a 95% MATA-Wald confidence interval, an' 
#' alpha=0.05 to a 90% interval.  'alpha' must be between 0 and 0.5.
#' Default value is alpha=0.025.
#' @param normal.lm Specify normal.lm=TRUE for the normal linear model case, and 
#' normal.lm=FALSE otherwise.  When normal.lm=TRUE, the argument 
#' 'residual.dfs' must also be supplied.  See USAGE section, 
#' and Turek and Fletcher (2012) for additional details.
#' @param residual.dfs A vector containing the residual (error) degrees of freedom 
#' under each candidate model.  This argument must be provided 
#' when the argument normal.lm=TRUE.
#' @param ... additional arguments passed to specific functions
#' @return If vcv=FALSE, the return value is a dataframe of model averaged
#' estimates and standard errors for a particular type of real parameter (e.g.,
#' Phi).  The design data are appended to the dataframe to enable subsettting
#' of the estimates based on features of the design data such as age, time,
#' cohort and grouping variables.
#' 
#' If vcv=TRUE, confidence interval (lcl,ucl) limits are added to the dataframe
#' which is contained in a list with the var-cov matrix.
#' @author Jeff Laake
#' @seealso \code{\link{collect.models}}, \code{\link{covariate.predictions}},
#' \code{\link{model.table}}, \code{\link{compute.links.from.reals}},
#' \code{\link{model.average.list}}
#' @references Burnham, K. P. and D. R. Anderson. 2002. Model Selection and
#' Multimodel Inference: A Practical Information-Theoretic Approach, Second
#' edition. Springer, New York.
#' @keywords utility
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(dipper)
#' run.dipper=function()
#' {
#' #
#' # Process data
#' #
#' dipper.processed=process.data(dipper,groups=("sex"))
#' #
#' # Create default design data
#' #
#' dipper.ddl=make.design.data(dipper.processed)
#' #
#' # Add Flood covariates for Phi and p that have different values
#' #
#' dipper.ddl$Phi$Flood=0
#' dipper.ddl$Phi$Flood[dipper.ddl$Phi$time==2 | dipper.ddl$Phi$time==3]=1
#' dipper.ddl$p$Flood=0
#' dipper.ddl$p$Flood[dipper.ddl$p$time==3]=1
#' #
#' #  Define range of models for Phi
#' #
#' Phi.dot=list(formula=~1)
#' Phi.time=list(formula=~time)
#' Phi.sex=list(formula=~sex)
#' Phi.sextime=list(formula=~sex+time)
#' Phi.sex.time=list(formula=~sex*time)
#' Phi.Flood=list(formula=~Flood)
#' #
#' #  Define range of models for p
#' #
#' p.dot=list(formula=~1)
#' p.time=list(formula=~time)
#' p.sex=list(formula=~sex)
#' p.sextime=list(formula=~sex+time)
#' p.sex.time=list(formula=~sex*time)
#' p.Flood=list(formula=~Flood)
#' #
#' # Collect pairings of models
#' #
#' cml=create.model.list("CJS")
#' #
#' # Run and return the list of models
#' #
#' return(mark.wrapper(cml,data=dipper.processed,ddl=dipper.ddl))
#' }
#' dipper.results=run.dipper()
#' Phi.estimates=model.average(dipper.results,"Phi",vcv=TRUE)
#' p.estimates=model.average(dipper.results,"p",vcv=TRUE)
#' run.dipper=function()
#' {
#' data(dipper)
#' dipper$nsex=as.numeric(dipper$sex)-1
#' #NOTE:  This generates random valules for the weights so the answers using
#' #  ~weight will vary
#' dipper$weight=rnorm(294)  
#' mod1=mark(dipper,groups="sex",
#'   model.parameters=list(Phi=list(formula=~sex+weight)))
#' mod2=mark(dipper,groups="sex",
#'   model.parameters=list(Phi=list(formula=~sex)))
#' mod3=mark(dipper,groups="sex",
#'   model.parameters=list(Phi=list(formula=~weight)))
#' mod4=mark(dipper,groups="sex",
#'   model.parameters=list(Phi=list(formula=~1)))
#' dipper.list=collect.models()
#' return(dipper.list)
#' }
#' dipper.results=run.dipper()
#' real.averages=model.average(dipper.results,vcv=TRUE)
#' # get model averaged estimates for all parameters and use average 
#' # covariate values in models with covariates
#' real.averages$estimates 
#' # get model averaged estimates for Phi using a value of 2 for weight
#' model.average(dipper.results,"Phi", 
#'   data=data.frame(weight=2),vcv=FALSE)  
#' # what you can't do yet is use different covariate values for 
#' # different groups to get covariances of estimates based on different
#' # covariate values; for example, you can get average survival of females 
#' # at average female weight and average survival of males at average 
#' # male weight in separate calls to model.average but not in the same call 
#' # to get covariances; however, if you standardized weight by group 
#' # (ie stdwt = weight - groupmean) then using 0 for the covariate value would give
#' # the model averaged Phi by group at the average group weights and its 
#' # covariance. You can do the above for
#' # a single model with find.covariates/fill.covariates.
#' # get model averaged estimates of first Phi(1) and first p(43) and v-c matrix
#' model.average(dipper.results,vcv=TRUE,indices=c(1,43))  
#' }
#' 
model.average.marklist<- function(x,parameter=NULL,data=NULL,vcv=FALSE,drop=TRUE,indices=NULL,revised=TRUE,mata=FALSE, normal.lm=FALSE, residual.dfs=0, alpha=0.025,...)
{
#
# Check validity of arguments and set up some variables
#
model.list=x
model=load.model(model.list[[1]])
if(class(model.list)!="marklist")
  stop("\nArgument for model.average must be a marklist created by collect.models\n")
if(is.null(model.list$model.table))
   stop("\nmarklist created by collect.models must contain a model.table to use model.average\n")
if(!is.null(parameter))
{
   if(!is.null(indices))
      message("\nNote: indices value has been ignored because parameter was set\n")
   if(!parameter%in%names(model$parameters))
      stop(paste("\n",parameter,"is not a valid parameter for these results\n"))
   parameters=setup.parameters(model$model)
   parameter.names=names(parameters)
   type=parameters[[match(parameter,parameter.names)]]$type
   begin.index=model$pims[[parameter]][[1]]$pim[1,1]
}
else
{
   begin.index=model$pims[[1]][[1]]$pim[1,1]
}
model.table=model.list$model.table
reals=vector("list",length=dim(model.table)[1])
#
# Determine if any of the models should be dropped because beta.var non-positive
#
dropped.models=NULL
for (i in 1:dim(model.table)[1])
{
  model=load.model(model.list[[i]])
  if(is.null(model$output) || is.null(model$results))
  {
     message("\nModel ",i,"is missing results and was excluded from model averaging\n")
     dropped.models=c(dropped.models,i)
  }
}
if(drop)
{
   for (i in 1:dim(model.table)[1])
   {
      model=load.model(model.list[[i]])
      if(!is.null(parameter))
      {
         indices=NULL
         for (j in 1:length(model$pims[[parameter]]))
           indices=c(indices,as.vector(t(model$pims[[parameter]][[j]]$pim)))
         indices=indices[indices>0]
      }
      model.indices=unique(model$simplify$pim.translation[indices])
      used.beta=which(apply(model$design.matrix[model.indices,,drop=FALSE],2,function(x)!all(x=="0")))
      if(any(diag(model$results$beta.vcv[used.beta,used.beta,drop=FALSE])<0))
#      if(any(diag(model$results$beta.vcv)<0))
      {
         dropped.models=c(dropped.models,i)
         message("\nModel ",i,"dropped from the model averaging because one or more beta variances are not positive\n")
      }
   }
#
# If any models have been dropped, recompute the weights for the model table
#
   if(!is.null(dropped.models))
   {
      model.table$weight[as.numeric(row.names(model.table))%in%dropped.models]=0
      model.table$weight=model.table$weight/sum(model.table$weight)
   }
}
#
#  Loop over each model and compute model averaged estimates and
#  model averaged correlation matrix
#
firstmodel=TRUE
for (i in 1:dim(model.table)[1])
{
   if(i %in% dropped.models) next
   model=load.model(model.list[[i]])
   if(is.null(parameter))
   {
      if(is.null(data))
      {
         if(!is.null(model$results$covariate.values$Value))
         {
            Covdata=as.data.frame(matrix(model$results$covariate.values$Value,nrow=1))
         } else
         {
            Covdata=NULL
         }
         if(!is.null(Covdata))
         {
            names(Covdata)=model$covariates
            zlist=compute.real(model,data=Covdata,se=TRUE,vcv=vcv)
         }
         else
            zlist=compute.real(model,design=model$design.matrix,se=TRUE,vcv=vcv)
      }
      else
         zlist=compute.real(model,data=data,se=TRUE,vcv=vcv)
      if(vcv)
      {
         z=list(estimates=data.frame(estimate=zlist$real,se=zlist$se.real,lcl=zlist$lcl,ucl=zlist$ucl,fixed=zlist$fixed),
                     vcv.real=zlist$vcv.real)
         z$estimates$par.index=1:dim(z$vcv.real)[1]
         row.names(z$vcv.real)=z$estimates$par.index
         colnames(z$vcv.real)=z$estimates$par.index
         if(!is.null(model$simplify))
         {
           z$estimates=z$estimates[model$simplify$pim.translation,]
           rownames(z$estimates)=model$simplify$real.labels
         }
      } else
      {
         z=zlist
         if(!is.null(model$simplify))
         {
           z=z[model$simplify$pim.translation,]
           rownames(z)=model$simplify$real.labels
         }
      }
      if(!is.null(indices))
      {
         if(vcv)
            z$estimates=z$estimates[indices,]
         else
            z=z[indices,]
      }
   }
   else
   {
      if(is.null(data))
      {
         if(!is.null(model$results$covariate.values$Value))
         {
            Covdata=as.data.frame(matrix(model$results$covariate.values$Value,nrow=1))
         } else
         {
            Covdata=NULL
         }
         if(!is.null(Covdata))
         {
            names(Covdata)=model$covariates
            z=get.real(model,parameter,data=Covdata,se=TRUE,vcv=vcv)
         }
         else
            z=get.real(model,parameter,design=model$design.matrix,se=TRUE,vcv=vcv)
      }
      else
         z=get.real(model,parameter,data=data,se=TRUE,vcv=vcv)
   }
   if(vcv)
      reals[[i]]=z$estimates
   else
      reals[[i]]=z
   if(firstmodel)
   {
      firstmodel=FALSE
      nreals=dim(reals[[i]])[1]
      estimates.average=rep(0,nreals)
      if(vcv) cor.average=matrix(0,nrow=nreals,ncol=nreals)
   }
   else
      if(dim(reals[[i]])[1]!=length(estimates.average))
         stop("\nCannot model average models with different structures\n")
   estimates.average=estimates.average+reals[[i]]$estimate*model.table$weight[as.numeric(row.names(model.table))==i]
   if(vcv)
   {
      rn=as.numeric(row.names(z$vcv.real))
      expanded.vcv=matrix(NA,nrow=max(rn),ncol=max(rn))
      zz=matrix(1,nrow=length(rn),ncol=length(rn))*rn
      expanded.vcv[cbind(as.vector(zz),as.vector(t(zz)))]=z$vcv.real
      xx=matrix(1,nrow=length(z$estimates$par.index),ncol=length(z$estimates$par.index))*z$estimates$par.index
      full.vcv=matrix(expanded.vcv[cbind(as.vector(xx),as.vector(t(xx)))],nrow=nreals,ncol=nreals)/
             outer(reals[[i]]$se,reals[[i]]$se,"*")
      full.vcv[is.nan(full.vcv)]=0
      full.vcv[is.infinite(full.vcv)]=0
	  if(any(is.infinite(diag(full.vcv)))) 
		  warning("Infinite correlation (se=0) for model  ",i, " for estimate ",which(is.infinite(diag(full.vcv))),"\n")
	  diag(full.vcv)=1
      cor.average=cor.average+full.vcv*model.table$weight[as.numeric(row.names(model.table))==i]
   }
}
#
#  Next compute model averaged se and create dataframe of estimates
#
se.average=rep(0,nreals)
for (i in 1:dim(model.table)[1])
{
   if(i %in% dropped.models) next
   if(revised)
      se.average=se.average+model.table$weight[as.numeric(row.names(model.table))==i]*
              (reals[[i]]$se^2 + (reals[[i]]$estimate-estimates.average)^2)  
   else
      se.average=se.average+model.table$weight[as.numeric(row.names(model.table))==i]*
              sqrt(reals[[i]]$se^2 + (reals[[i]]$estimate-estimates.average)^2)
}
if(revised) se.average=sqrt(se.average)
first.index=((1:dim(model.table)[1])[!(1:dim(model.table)[1])%in%dropped.models])[1]
if(!is.null(parameter))
   other.values=summary(model.list[[first.index]],se=TRUE)$reals[[parameter]]
else
   other.values=summary(model.list[[first.index]],se=TRUE)$reals

if(!vcv)
{
   if(is.null(parameter))
      if(is.null(indices))
         result=cbind(data.frame(par.index=begin.index:(begin.index+nreals-1),estimate=estimates.average,se=se.average))
      else
         result=cbind(data.frame(par.index=indices,estimate=estimates.average,se=se.average))
   else
      result=cbind(data.frame(par.index=begin.index:(begin.index+nreals-1),estimate=estimates.average,se=se.average),other.values[,(7:dim(other.values)[2])])
   return(result)
}
#
# If vcv requested then compute model averaged vc from model averaged correlations
# and model averaged conf intervals for the parameters
#
else
{
   if(is.null(parameter)&!is.null(indices))
      result=cbind(data.frame(par.index=indices,estimate=estimates.average,se=se.average))
   else
      result=cbind(data.frame(par.index=begin.index:(begin.index+nreals-1),estimate=estimates.average,se=se.average))
   vcv.real=cor.average*outer(se.average,se.average,"*")
   vcv.real[is.nan(vcv.real)]=0
   vcv.real[is.infinite(abs(vcv.real))]=0  
   row.names(vcv.real)=result$par.index
   colnames(vcv.real)=row.names(vcv.real)
   link.list=compute.links.from.reals(result$estimate,model.list[[1]],parm.indices=result$par.index,vcv.real=vcv.real,use.mlogits=FALSE)
   if(!mata)
   {
	   result$lcl=link.list$estimates-qnorm(1-alpha,0,1)*sqrt(diag(link.list$vcv))
	   result$ucl=link.list$estimates+qnorm(1-alpha,0,1)*sqrt(diag(link.list$vcv))
	   result$lcl=apply(data.frame(x=result$lcl,links=link.list$links),1,function(x){inverse.link(as.numeric(x[1]),x[2])})
	   result$ucl=apply(data.frame(x=result$ucl,links=link.list$links),1,function(x){inverse.link(as.numeric(x[1]),x[2])})
	   result$lcl[is.na(result$lcl)]=result$estimate[is.na(result$lcl)]
	   result$ucl[is.na(result$ucl)]=result$estimate[is.na(result$ucl)]
   } else
   {
	   weights=NULL
	   for (j in 1:length(reals))
		   if(!is.null(reals[[j]]))weights=c(weights,model.table$weight[j])
	   for(i in 1:nrow(result))
	   {
		   link.estimate=NULL
		   link.se=NULL
		   for (j in 1:length(reals))
		   {
			   if(!is.null(reals[[j]]))
			   {
				   link.list=compute.links.from.reals(reals[[j]]$estimate[i],model.list[[1]],parm.indices=result$par.index[i],vcv.real=vcv.real[result$par.index[i],result$par.index[i]],use.mlogits=FALSE)
				   link.estimate=c(link.estimate,link.list$estimates)
				   link.se=c(link.se,sqrt(diag(link.list$vcv)))
			   }	   
		   }
		   interval=mata.wald(theta.hats=link.estimate, se.theta.hats=link.se, model.weights=weights, normal.lm=normal.lm, residual.dfs=residual.dfs, alpha=alpha) 	 
		   result$lcl[i]=inverse.link(interval[1],link.list$links)
		   result$ucl[i]=inverse.link(interval[2],link.list$links)
	   }
	   result$lcl[is.na(result$lcl)]=result$estimate[is.na(result$lcl)]
	   result$ucl[is.na(result$ucl)]=result$estimate[is.na(result$ucl)]
   }
   if(!is.null(parameter)) result=cbind(result,other.values[,(7:dim(other.values)[2])])
   return(list(estimates=result,vcv.real=vcv.real))
}
}







