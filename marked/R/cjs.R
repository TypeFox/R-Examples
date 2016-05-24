#' Fitting function for CJS models
#' 
#' A function for computing MLEs for a specified Cormack-Jolly-Seber open
#' population capture-recapture model for processed dataframe \code{x} with
#' user specified formulas in \code{parameters} that create list of design
#' matrices \code{dml}. This function can be called directly but is most easily
#' called from \code{\link{crm}} that sets up needed arguments.
#' 
#' It is easiest to call \code{cjs} through the function \code{\link{crm}}.
#' Details are explained there.
#' 
#' Be cautious with this function at present.  It does not include many checks
#' to make sure values like fixed values will remain in the specified range of
#' the data.  Normally this would not be a big problem but because
#' \code{\link{cjs.lnl}} calls an external FORTRAN subroutine, if it gets a
#' subscript out of bounds, it will cause R to terminate.  So make sure to save
#' your workspace frequently if you use this function in its current
#' implementation.
#' 
#' @param x processed dataframe created by process.data
#' @param ddl list of dataframes for design data; created by call to
#' \code{\link{make.design.data}}
#' @param dml list of design matrices created by \code{\link{create.dm}} from
#' formula and design data
#' @param model_data a list of all the relevant data for fitting the model including
#' imat, Phi.dm,p.dm,Phi.fixed,p.fixed, and time.intervals. It is used to save values
#' and avoid accumulation again if the model was re-rerun with an additional call to cjs when
#' using autoscale or re-starting with initial values.  It is stored with returned model object.
#' @param parameters equivalent to \code{model.parameters} in \code{\link{crm}}
#' @param accumulate if TRUE will accumulate capture histories with common
#' value and with a common design matrix for Phi and p to speed up execution
#' @param initial list of initial values for parameters if desired; if each is a named vector
#' from previous run it will match to columns with same name
#' @param method method to use for optimization; see \code{optim}
#' @param hessian if TRUE will compute and return the hessian
#' @param debug if TRUE will print out information for each iteration
#' @param chunk_size specifies amount of memory to use in accumulating capture
#' histories; amount used is 8*chunk_size/1e6 MB (default 80MB)
#' @param refit non-zero entry to refit
#' @param itnmax maximum number of iterations
#' @param control control string for optimization functions
#' @param scale vector of scale values for parameters
#' @param use.admb if TRUE creates data file for admbcjs.tpl and runs admb optimizer
#' @param crossed if TRUE it uses cjs.tpl or cjs_reml.tpl if reml=FALSE or TRUE respectively; if FALSE, then it uses cjsre which can use Gauss-Hermite integration
#' @param compile if TRUE forces re-compilation of tpl file
#' @param extra.args optional character string that is passed to admb if use.admb==TRUE
#' @param reml if set to TRUE uses cjs_reml if crossed 
#' @param clean if TRUE, deletes the tpl and executable files for amdb if use.admb=T
#' @param ... any remaining arguments are passed to additional parameters
#' passed to \code{optim} or \code{\link{cjs.lnl}}
#' @import R2admb optimx
#' @return The resulting value of the function is a list with the class of
#' crm,cjs such that the generic functions print and coef can be used.
#' Elements are 1) beta: named vector of parameter estimatesm 2) lnl: -2*log
#' likelihood, 3) AIC: lnl + 2* number of parameters, 4) convergence: result from \code{optim}; if 0 \code{optim} thinks it
#' converged, 5) count:\code{optim} results of number of function
#' evaluations, 6) reals: dataframe of data and real Phi and p estimates for
#' each animal-occasion excluding those that occurred before release, 7) vcv:var-cov matrix of betas if hessian=TRUE was set.
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @references Pledger, S., K. H. Pollock, et al. (2003). Open
#' capture-recapture models with heterogeneity: I. Cormack-Jolly-Seber model.
#' Biometrics 59(4):786-794.
cjs=function(x,ddl,dml,model_data=NULL,parameters,accumulate=TRUE,initial=NULL,method,
            hessian=FALSE,debug=FALSE,chunk_size=1e7,refit,itnmax=NULL,control=NULL,scale,
			use.admb=FALSE,crossed=TRUE,compile=FALSE,extra.args=NULL,reml,clean=TRUE,...)
{
   if(use.admb)accumulate=FALSE
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
#  Store data from x into x
   x=x$data
#  set default frequencies if not used
   freq=NULL
   if(!is.null(x$freq))freq=x$freq
#  get first and last vectors, loc and chmat with process.ch and store in imat
   ch=x$ch
   imat=process.ch(ch,freq,all=FALSE)
#  Use specified initial values or create if null
   if(is.null(initial))
	   par=cjs.initial(dml,imat)
   else
       par=set.initial(names(dml),dml,initial)$par
   initial=par
#  Create list of model data for optimization
	model_data=list(Phi.dm=dml$Phi$fe,p.dm=dml$p$fe,imat=imat,Phi.fixed=parameters$Phi$fixed,
			p.fixed=parameters$p$fixed,time.intervals=time.intervals)
#   If data are to be accumulated based on ch and design matrices do so here;
#   Problems with accumulation and fixed values 10 Jan 2014; turned off accumulate if fixed
    if(parameters$p$fixed[1,1]>0 | parameters$Phi$fixed[1,1]>0) accumulate=FALSE
	if(accumulate)
	{
		message("Accumulating capture histories based on design. This can take awhile...")
		#flush.console()
		model_data.save=model_data   
		model_data=cjs.accumulate(x,model_data,nocc,freq,chunk_size=chunk_size)
	}else
		model_data.save=model_data
#   Create links  -- not used at present; idea here is to use sin links for parameters where you can   
#   Phi.links=create.links(Phi.dm)
#   Phi.links=which(Phi.links==1)
#   p.links=create.links(p.dm)
#   p.links=which(p.links==1)
#  Scale the design matrices and parameters with either input scale or computed scale
   if(use.admb)scale=1
   scale=set.scale(names(dml),model_data,scale)
   model_data=scale.dm(model_data,scale)
########################################################################
#   CJS with R
########################################################################
   if(!use.admb)
   {
	   par=scale.par(par,scale)
	   #  Call optimx to find mles with cjs.lnl which gives -log-likelihood
	   message("Starting optimization for ",length(par)," parameters...")
	   #flush.console()
	   markedfunc_eval=0
	   cjsenv=environment()
	   if("SANN"%in%method)
	   {
		   mod=optim(par,cjs.lnl,model_data=model_data,method="SANN",hessian=FALSE,
				   debug=debug,control=control,cjsenv=cjsenv,...)
		   par= mod$par
		   convergence=mod$convergence
		   lnl=mod$value
		   counts=mod$counts
	   }else
	   {
		   mod=optimx(par,cjs.lnl,model_data=model_data,method=method,hessian=FALSE,
						   debug=debug,control=control,itnmax=itnmax,cjsenv=cjsenv,...)
		   par <- coef(mod, order="value")[1, ]
		   mod=as.list(summary(mod, order="value")[1, ])
		   convergence=mod$convcode
		   lnl=mod$value
	   }
	   #  Rescale parameter vector 
	   cjs.beta=unscale.par(par,scale)
       # Create results list 
	   res=list(beta=cjs.beta,neg2lnl=2*lnl,AIC=2*lnl+2*sum(sapply(cjs.beta,length)),
			   convergence=convergence,optim.details=mod,
			   model_data=model_data,
			   options=list(scale=scale,accumulate=accumulate,initial=initial,method=method,
					   chunk_size=chunk_size,itnmax=itnmax,control=control))
       # Compute hessian if requested
	   if(hessian) 
	   {
		   message("Computing hessian...")
		   res$beta.vcv=cjs.hessian(res)
	   } 
   } else
   {
########################################################################
#      CJS with ADMB
########################################################################
	   if(R.Version()$os=="mingw32")
		   ext=".exe"
       else
		   ext=""
	   sdir=system.file(package="marked")
	   # set tpl filename
	   if(!crossed)
	   {
		   tpl="cjsre"
	   } else
	   {
		   if(reml) 
			   tpl="cjs_reml"
		   else
			   tpl="cjs"
	   }
	   # cleanup any leftover files
	   clean_admb(tpl)
	   # if argument clean is TRUE, delete exe and TPL files as well
	   if(clean)
	   {
		  if(file.exists(paste(tpl,".tpl",sep=""))) unlink(paste(tpl,".tpl",sep=""))
		  if(file.exists(paste(tpl,ext,sep=""))) unlink(paste(tpl,ext,sep=""))
	  }
	   # if tpl is not available, copy from the package directory
	   if(!file.exists(paste(tpl,".tpl",sep="")))
	   {
		   file.copy(file.path(sdir,paste(tpl,".tpl",sep="")),file.path(getwd(),paste(tpl,".tpl",sep="")),overwrite=TRUE)
		   file.copy( file.path(sdir,list.files(sdir,"*.cpp")),file.path(getwd()),overwrite=TRUE)
	   } 
	   # if exe available in the package or workspace directory use it
       have.exe=TRUE
	   if(!file.exists(file.path(getwd(),paste(tpl,ext,sep=""))))
	   {
		   if(file.exists(file.path(sdir,paste(tpl,ext,sep=""))) &!compile)
		   {
			   file.copy(file.path(sdir,paste(tpl,ext,sep="")),file.path(getwd(),paste(tpl,ext,sep="")),overwrite=TRUE)
		   }else
 	       # if there is no exe in either place then check to make sure ADMB is installed and linked
		   {
			  have.exe=FALSE
		   }
	   }
	   if(!have.exe | compile)
	   {
		   # no exe or compile set TRUE; see if admb can be found; this is not a complete test but should catch the novice user who has
	       # not setup admb at all
	       if(R.Version()$os=="mingw32" & Sys.which(paste("tpl2cpp",ext,sep=""))=="")
			  stop("admb not found; setup links to admb and c++ compiler with environment variables or put in path")
		   else
		   {
			  compile_admb(tpl,re=TRUE,safe=TRUE,verbose=T)
		   }
	   }
	   # create admbcjs.dat file to create its contents 
	   con=file(paste(tpl,".dat",sep=""),open="wt")
	   # Number of observations
	   n=length(model_data$imat$freq)
	   write(n,con,append=FALSE)
	   # Number of occasions
	   nocc=model_data$imat$nocc
	   write(nocc,con,append=TRUE)
	   write(as.numeric(debug),con,append=TRUE)
	   # capture history matrix
	   write(t(model_data$imat$chmat),con,ncolumns=nocc,append=TRUE)
	   # first occasions seen 
	   write(model_data$imat$first,con,ncolumns=n,append=TRUE)
	   # last occasions seen 
	   write(model_data$imat$last,con,ncolumns=n,append=TRUE)
	   # frequency of capture history 
       if(tpl=="cjsre") write(model_data$imat$freq,con,ncolumns=n,append=TRUE)
	   # indicator for loss on capture 
	   write(model_data$imat$loc,con,ncolumns=n,append=TRUE)
	   write(t(model_data$time.intervals),con,ncolumns=nocc-1,append=TRUE)
       #phi dm portion
       phidm=model_data$Phi.dm
	   phidm=cbind(phidm,rep(-1,nrow(phidm)))
	   if(model_data$Phi.fixed[1,1]!= -1)
	   {
		   index=(nocc-1)*(model_data$Phi.fixed[,1]-1)+model_data$Phi.fixed[,2]
		   phidm[index,ncol(phidm)]=model_data$Phi.fixed[,3]
		   phidm[index,1:(ncol(phidm)-1)]=0
	   }
	   write(ncol(phidm)-1,con,append=TRUE)
	   write(t(as.matrix(phidm)),con,ncolumns=ncol(phidm),append=TRUE)
	   phimixed=mixed.model.admb(parameters$Phi$formula,ddl$Phi)
       nphisigma=0
	   if(!is.null(phimixed$re.dm))nphisigma=ncol(phimixed$re.dm)
	   if(crossed & !is.null(phimixed$re.dm))
       {
		   phimixed$re.indices[ddl$Phi$Time<ddl$Phi$Cohort,]=NA
		   phimixed=reindex(phimixed,ddl$Phi$id)
       }
	   mixed.model.dat(phimixed,con,!crossed,n)
	   # p dm portion
       pdm=model_data$p.dm
	   pdm=cbind(pdm,rep(-1,nrow(pdm)))
	   if(model_data$p.fixed[1,1]!= -1)
	   {
		   index=(nocc-1)*(model_data$p.fixed[,1]-1)+model_data$p.fixed[,2]-1
		   pdm[index,ncol(pdm)]=model_data$p.fixed[,3]
		   pdm[index,1:(ncol(pdm)-1)]=0
	   }
	   write(ncol(pdm)-1,con,append=TRUE)
       write(t(as.matrix(pdm)),con,ncolumns=ncol(pdm),append=TRUE)
       pmixed=mixed.model.admb(parameters$p$formula,ddl$p)
	   npsigma=0
	   if(!is.null(pmixed$re.dm))npsigma=ncol(pmixed$re.dm)
	   if(crossed &!is.null(pmixed$re.dm))
       {
		   pmixed$re.indices[ddl$p$Time<ddl$p$Cohort,]=NA
		   pmixed=reindex(pmixed,ddl$p$id)
       }	   
       mixed.model.dat(pmixed,con,!crossed,n)
	   close(con)
       # setup initial value file (.pin)
	   con=file(paste(tpl,".pin",sep=""),open="wt")
	   write(par$Phi,con,ncolumns=length(par$Phi),append=FALSE)
	   write(par$p,con,ncolumns=length(par$p),append=TRUE)
	   if(is.null(extra.args)) extra.args=""
	   if(nphisigma+npsigma>0) 
	   {
		   warning("\nReal parameter estimates are not produced currently for random effect models\n") 
		   if(is.null(initial$sigmaPhi))
			   sigma.initial=rep(-2,nphisigma)
		   else
		   {
			   sigma.initial=initial$sigmaPhi
			   if(length(initial$sigmaPhi)!=nphisigma)stop("length of initial values for sigmaPhi does not match design")
		   }
		   if(is.null(initial$sigmap))
			   sigma.initial=c(sigma.initial,rep(-2,npsigma))
		   else
		   {
			   sigma.initial=c(sigma.initial,initial$sigmap)
			   if(length(initial$sigmap)!=npsigma)stop("length of initial values for sigmap does not match design")
		   }
		   write(sigma.initial,con,ncolumns=1,append=TRUE)
	  	   if(nphisigma>0) write(rep(0,nphisigma*nrow(phimixed$re.dm)),con,ncolumns=nphisigma,append=TRUE)
		   if(npsigma>0) write(rep(0,npsigma*nrow(pmixed$re.dm)),con,ncolumns=npsigma,append=TRUE)
		   if(crossed)
		   {
			   if(any(model_data$imat$freq!=1)) stop("\n freq cannot be > 1 if crossed effects; don't accumulate")
			   extra.args=paste(extra.args,"-shess")
		   } else
			   extra.args=paste(extra.args,"-gh 14")
	   }
	   close(con)   
	   cat("\nrunning ADMB program\n")
	   flush.console()
	   if(!hessian)extra.args=paste(extra.args,"-nohess")
	   xx=run_admb(tpl,verbose=T,extra.args=extra.args)
	   convergence=attr(xx,"status")
	   if(is.null(convergence))convergence=0
	   res=read_admb(tpl)
	   cjs.beta.fixed=unscale.par(c(res$coeflist$phi_beta,res$coeflist$p_beta),scale)
	   cjs.beta.random=c(Phi=res$coeflist$phi_sigma,p=res$coeflist$p_sigma)
	   if(!is.null(cjs.beta.random))names(cjs.beta.random)=paste("sigma_",names(cjs.beta.random),sep="")
	   cjs.beta=c(cjs.beta.fixed,cjs.beta.random)
	   parnames=names(unlist(cjs.beta))
	   fixed.npar=length(cjs.beta.fixed)
       if(fixed.npar<res$npar)
	   {
		   allnames=names(unlist(res$coeflist))[1:(res$npar+res$npar_re)]
		   allnames=sub("phi_","Phi.",allnames)
		   allnames=sub("p_","p.",allnames)
		   allnames[1:fixed.npar]=parnames[1:fixed.npar]
		   random.effects=coef(res)[(fixed.npar+1):res$npar]
		   names(random.effects)=allnames[(fixed.npar+1):res$npar]
	   }else
	   {
		   random.effects=NULL
		   allnames=parnames
	   }
	   beta=list(cjs.beta)
	   if(!is.null(res$hes))
	   {
		   beta.vcv=solvecov(res$hes)$inv
		   rownames(res$hes)=allnames
		   colnames(res$hes)=rownames(res$hes)
		   if(all(diag(beta.vcv>0))) 
		      res$cor=beta.vcv/outer(sqrt(diag(beta.vcv)),sqrt(diag(beta.vcv)))
	   }  else
		   beta.vcv=res$vcov
	   rownames(beta.vcv)=allnames
	   colnames(beta.vcv)=rownames(beta.vcv)
	   rownames(res$cor)=rownames(beta.vcv)
	   colnames(res$cor)=rownames(beta.vcv)
	   res$vcov=NULL
	   optim.details=c(fn=res$fn,maxgrad=res$maxgrad,eratio=res$eratio)
	   options=list(extra.args=extra.args)
	   res$cor=NULL
	   res$maxgrad=NULL
	   results=c(beta=beta,neg2lnl=-2*res$loglik,AIC=-2*res$loglik+2*res$npar,convergence=convergence)
	   results$optim.details=optim.details
	   results$options=options
	   results$random.effects=res$random.effects
	   results$coeflist=res$coeflist
	   results$npar=list(npar=res$npar,npar_re=res$npar_re,npar_sdrpt=res$npar_sdrpt,npar_total=res$npar_total)
	   results$beta.vcv=beta.vcv
	   res=results
   }
#  Restore non-accumulated, non-scaled dm's etc
   res$model_data=model_data.save
#  Assign S3 class values and return
   class(res)=c("crm","mle","cjs")
   return(res)
}


