#' Fitting function for Jolly-Seber model using Schwarz-Arnason POPAN
#' formulation
#' 
#' A function for computing MLEs for a specified Jolly-Seber open population
#' capture-recapture model for processed dataframe \code{x} with user specified
#' formulas in \code{parameters} that create list of design matrices
#' \code{dml}. This function can be called directly but is most easily called
#' from \code{\link{crm}} that sets up needed arguments.
#' 
#' It is easiest to call \code{js} through the function \code{\link{crm}}.
#' Details are explained there.
#' 
#' Be cautious with this function at present.  It does not include many checks
#' to make sure values like fixed values will remain in the specified range of
#' the data.  Normally this would not be a big problem but because
#' \code{\link{js.lnl}} calls an external FORTRAN subroutine via
#' \code{\link{cjs.lnl}}, if it gets a subscirpt out of bounds, it will cause R
#' to terminate.  So make sure to save your workspace frequently if you use
#' this function in its current implementation.
#' 
#' @param x processed dataframe created by process.data
#' @param ddl list of dataframes for design data; created by call to
#' \code{\link{make.design.data}}
#' @param dml list of design matrices created by \code{\link{create.dm}} from
#' formula and design data
#' @param model_data a list of all the relevant data for fitting the model including
#' imat, Phi.dm,p.dm,Phi.fixed,p.fixed, and time.intervals. It is used to save values
#' and avoid accumulation again if the model was re-rerun with an additional call to js when
#' using autoscale or re-starting with initial values.  It is stored with returned model object.
#' @param parameters equivalent to \code{model.parameters} in \code{\link{crm}}
#' @param accumulate if TRUE will accumulate capture histories with common
#' value and with a common design matrix for Phi and p to speed up execution
#' @param initial initial values for parameters if desired; if named vector
#' from previous run it will match to columns with same name
#' @param method method to use for optimization; see \code{optimx}
#' @param hessian if TRUE will compute and return the hessian
#' @param debug if TRUE will print out information for each iteration
#' @param chunk_size specifies amount of memory to use in accumulating capture
#' histories; amount used is 8*chunk_size/1e6 MB (default 80MB)
#' @param refit non-zero entry to refit
#' @param itnmax maximum number of iterations
#' @param control control string for optimization functions
#' @param scale vector of scale values for parameters
#' @param ... any remaining arguments are passed to additional parameters
#' passed to \code{optimx} or \code{\link{js.lnl}}
#' @return The resulting value of the function is a list with the class of
#' crm,js such that the generic functions print and coef can be used.
#' \item{beta}{named vector of parameter estimates} \item{lnl}{-2*log
#' likelihood} \item{AIC}{lnl + 2* number of parameters}
#' \item{convergence}{result from \code{optimx}; if 0 \code{optimx} thinks it
#' converged} \item{count}{\code{optimx} results of number of function
#' evaluations} \item{reals}{dataframe of data and real Phi and p estimates for
#' each animal-occasion excluding those that occurred before release}
#' \item{vcv}{var-cov matrix of betas if hessian=TRUE was set}
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @references Schwarz, C. J., and A. N. Arnason. 1996. A general methodology
#' for the analysis of capture-recapture experiments in open populations.
#' Biometrics 52:860-873.
js=function(x,ddl,dml,model_data=NULL,parameters,accumulate=TRUE,initial=NULL,method="BFGS",
            hessian=FALSE,debug=FALSE,chunk_size=1e7,refit,itnmax=NULL,control=NULL,scale,...)
{
    nocc=x$nocc
#  Time intervals has been changed to a matrix (columns=intervals,rows=animals)
#  so that the initial time interval can vary by animal; use x$intervals if none are in ddl$Phi
	if(!is.null(ddl$Phi$time.interval))
		time.intervals=matrix(ddl$Phi$time.interval,nrow(x$data),ncol=nocc-1,byrow=TRUE)
	else
	if(is.vector(x$time.intervals))
		time.intervals=matrix(x$time.intervals,nrow=nrow(x$data),ncol=nocc-1,byrow=TRUE)
	else
		time.intervals=x$time.intervals
#  Create fixed matrices in parameters
   parameters=create.fixed.matrix(ddl,parameters)
   parameters$N$fixed=NULL
# Compute nobstot - total unique critters
   if(nrow(dml$N$fe)==1) 
      nobstot=sum(x$data$freq)    
   else
     nobstot=tapply(x$data$freq,x$data$group,sum)   
#  Store data from x into x
   x=x$data
#  set default frequencies if not used
   freq=NULL
   if(!is.null(x$freq))freq=x$freq
 # get first and last vectors, loc and chmat
   ch=x$ch
   imat=process.ch(ch,freq,all=TRUE)
#  Use specified initial values or create if null
	if(is.null(initial))
		initial=cjs.initial(dml,imat)
	par=set.initial(names(dml),dml,initial)$par
	initial=par
#  Create list of model data for optimization; if passed as an argument create model_data.save 
#  and use model_data (accumulated values); otherwise create model_data, save it and accumulate it
#  if requested.
   model_data=list(Phi.dm=dml$Phi$fe,p.dm=dml$p$fe,pent.dm=dml$pent$fe,N.dm=dml$N$fe,imat=imat,Phi.fixed=parameters$Phi$fixed,
		   p.fixed=parameters$p$fixed,pent.fixed=parameters$pent$fixed,time.intervals=time.intervals)
#  If data are to be accumulated based on ch and design matrices do so here;
   if(accumulate)
   {
 	   cat("Accumulating capture histories based on design. This can take awhile.\n")
	   flush.console()
	   model_data.save=model_data   
	   model_data=js.accumulate(x,model_data,nocc,freq,chunk_size=chunk_size)
   }else
	   model_data.save=NULL
#  Scale the design matrices and parameters with either input scale or computed scale
   scale=set.scale(names(dml),model_data,scale)
   model_data=scale.dm(model_data,scale)
   par=scale.par(par,scale)
#  call optim to find mles with js.lnl which gives -log-likelihood
   markedfunc_eval=0
   jsenv=environment()
   cat("Starting optimization",length(par)," parameters\n")
   if("SANN"%in%method)
   {
	   mod=optim(par,js.lnl,model_data=model_data,Phi.links=NULL,p.links=NULL,method="SANN",hessian=FALSE,
			   debug=debug,control=control,jsenv=jsenv,...)
	   par= mod$par
	   convergence=mod$convergence
	   lnl=mod$value
	   counts=mod$counts
   }else
   {
	   mod=optimx(par,js.lnl,model_data=model_data,method=method,hessian=hessian,
					   debug=debug,control=control,itnmax=itnmax,nobstot=nobstot,jsenv=jsenv,...)
	   par <- coef(mod, order="value")[1, ]
	   mod=as.list(summary(mod, order="value")[1, ])
	   convergence=mod$convcode
	   lnl=mod$value
   }
   js.beta=unscale.par(par,scale)
#  Compute additional likelihood component so lnl matches output for MARK; not needed for optimization
   if(is.null(model_data.save))
   {
	   if(is.null(x$group))
		   ui=tapply(model_data$imat$freq,list(model_data$imat$first),sum)
	   else
		   ui=tapply(model_data$imat$freq,list(model_data$imat$first,x$group),sum)	  
   } else
   {
	   if(is.null(x$group))
		   ui=tapply(model_data.save$imat$freq,list(model_data.save$imat$first),sum)
	   else
		   ui=tapply(model_data.save$imat$freq,list(model_data.save$imat$first,x$group),sum)
   }
   lnl=lnl+sum(lfactorial(ui))
#  create results list
   res=list(beta=js.beta,neg2lnl=2*lnl,AIC=2*lnl+2*sum(sapply(js.beta,length)),
		   convergence=convergence,optim.details=mod,
		   model_data=model_data,ns=nobstot,
		   options=list(scale=scale,accumulate=accumulate,initial=initial,method=method,
		   chunk_size=chunk_size,itnmax=itnmax,control=control))
#  Compute hessian if specified   
   if(hessian) 
   {
	   cat("Computing hessian\n")
	   res$beta.vcv=js.hessian(res)
   }   
#  Restore complete non-accumulated model_data with unscaled design matrices 
   if(!is.null(model_data.save)) res$model_data=model_data.save
#  Assign class and return
   class(res)=c("crm","mle","js")
   return(res)
}


