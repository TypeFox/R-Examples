#' Fitting function for Multistate CJS models
#' 
#' A function for computing MLEs for a Multi-state Cormack-Jolly-Seber open
#' population capture-recapture model for processed dataframe \code{x} with
#' user specified formulas in \code{parameters} that create list of design
#' matrices \code{dml}. This function can be called directly but is most easily
#' called from \code{\link{crm}} that sets up needed arguments.
#' 
#' It is easiest to call \code{mscjs} through the function \code{\link{crm}}.
#' Details are explained there.
#' 
#' @param x processed dataframe created by process.data
#' @param ddl list of dataframes for design data; created by call to
#' \code{\link{make.design.data}}
#' @param dml list of design matrices created by \code{\link{create.dm}} from
#' formula and design data
#' @param model_data a list of all the relevant data for fitting the model including
#' imat, S.dm,p.dm,Psi.dm,S.fixed,p.fixed,Psi.fixed and time.intervals. It is used to save values
#' and avoid accumulation again if the model was re-rerun with an additional call to cjs when
#' using autoscale or re-starting with initial values.  It is stored with returned model object.
#' @param parameters equivalent to \code{model.parameters} in \code{\link{crm}}
#' @param accumulate if TRUE will accumulate capture histories with common
#' value and with a common design matrix for S and p to speed up execution
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
#' @param re if TRUE creates random effect model admbcjsre.tpl and runs admb optimizer
#' @param compile if TRUE forces re-compilation of tpl file
#' @param extra.args optional character string that is passed to admb if use.admb==TRUE
#' @param clean if TRUE, deletes the tpl and executable files for amdb if use.admb=T
#' @param ... any remaining arguments are passed to additional parameters
#' passed to \code{optim} or \code{\link{cjs.lnl}}
#' @import R2admb
#' @return The resulting value of the function is a list with the class of
#' crm,cjs such that the generic functions print and coef can be used.
#' \item{beta}{named vector of parameter estimates} \item{lnl}{-2*log
#' likelihood} \item{AIC}{lnl + 2* number of parameters}
#' \item{convergence}{result from \code{optim}; if 0 \code{optim} thinks it
#' converged} \item{count}{\code{optim} results of number of function
#' evaluations} \item{reals}{dataframe of data and real S and p estimates for
#' each animal-occasion excluding those that occurred before release}
#' \item{vcv}{var-cov matrix of betas if hessian=TRUE was set}
#' @author Jeff Laake <jeff.laake@@noaa.gov>
#' @references Ford, J. H., M. V. Bravington, and J. Robbins. 2012. Incorporating individual variability into mark-recapture models. Methods in Ecology and Evolution 3:1047-1054.
#' 
mscjs=function(x,ddl,dml,model_data=NULL,parameters,accumulate=TRUE,initial=NULL,method,
		hessian=FALSE,debug=FALSE,chunk_size=1e7,refit,itnmax=NULL,control=NULL,scale,
		use.admb=FALSE,re=FALSE,compile=FALSE,extra.args="",clean=TRUE,...)
{
	accumulate=FALSE
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
#  Store data from x$data into x
	strata.labels=x$strata.labels
	uS=x$unobserved
	x=x$data
#  set default frequencies if not used
	freq=NULL
	if(!is.null(x$freq))freq=x$freq
#  get first and last vectors, loc and chmat with process.ch and store in imat
	ch=x$ch
	imat=process.ch(ch,freq,all=FALSE)
	chmat=matrix((unlist(strsplit(ch,","))),byrow=TRUE,ncol=nocc,nrow=length(ch))
	for(nlabel in 1:length(strata.labels))
		chmat=t(apply(chmat,1,sub,pattern=strata.labels[nlabel],replacement=nlabel))
#  Use specified initial values or create if null
	if(is.null(initial))
		par=list(Psi=rep(0,ncol(dml$Psi$fe)),
		         p=rep(0,ncol(dml$p$fe)),
		         S=rep(0,ncol(dml$S$fe)))
	else
		par=set.initial(names(dml),dml,initial)$par
#  Create list of model data for optimization
	model_data=list(S.dm=dml$S$fe,p.dm=dml$p$fe,Psi.dm=dml$Psi$fe,imat=imat,S.fixed=parameters$S$fixed,
			p.fixed=parameters$p$fixed,Psi.fixed=parameters$Psi$fixed,time.intervals=time.intervals)
#   If data are to be accumulated based on ch and design matrices do so here;
	if(accumulate)
	{
		cat("Accumulating capture histories based on design. This can take awhile.\n")
		flush.console()
		model_data.save=model_data   
		#model_data=mscjs.accumulate(x,model_data,nocc,freq,chunk_size=chunk_size)
	}else
		model_data.save=model_data
#   Create links  -- not used at present; idea here is to use sin links for parameters where you can   
#   S.links=create.links(S.dm)
#   S.links=which(S.links==1)
#   p.links=create.links(p.dm)
#   p.links=which(p.links==1)
#  Scale the design matrices and parameters with either input scale or computed scale
    use.admb=TRUE
    scale=1
	scale=set.scale(names(dml),model_data,scale)
	model_data=scale.dm(model_data,scale)
	if(R.Version()$os=="mingw32")
		ext=".exe"
	else
		ext=""
	# see if admb can be found; this is not a complete test but should catch the novice user who has
	# not setup admb at all
    if(Sys.which(paste("tpl2cpp",ext,sep=""))=="") 
	    stop("admb not found; setup links to admb and c++ compiler with environment variables or put in path")
	sdir=system.file(package="marked")
	# if multistate.tpl is not available copy from the package directory
	if(!re)
		tpl="multistate"
	else
	{
		tpl="notdone"
		if(!hessian)message("ignoring hessian setting; set to TRUE")
		hessian=TRUE
	}
    # cleanup any leftover admbcjs files
	clean_admb(tpl)
	# if argument clean is TRUE, delete exe and TPL files as well
	if(clean)
	{
		if(file.exists(paste(tpl,".tpl",sep=""))) unlink(paste(tpl,".tpl",sep=""))
		if(file.exists(paste(tpl,ext,sep=""))) unlink(paste(tpl,ext,sep=""))
	}
	if(!file.exists(paste(tpl,".tpl",sep="")))
		file.copy(file.path(sdir,paste(tpl,".tpl",sep="")),file.path(getwd(),paste(tpl,".tpl",sep="")),overwrite=TRUE)
	if(!file.exists(paste(tpl,ext,sep=""))|compile)
		compile_admb(tpl,re=re)
	# create .dat file to write its contents 
	con=file(paste(tpl,".dat",sep=""),open="wt")
	# Number of observations
	n=length(model_data$imat$freq)
	write(n,con,append=FALSE)
	# Number of occasions
	nocc=model_data$imat$nocc
	write(nocc,con,append=TRUE)
	# Number of states
	write(length(strata.labels),con,append=TRUE)
	# Number of unobserved states
	write(uS,con,append=TRUE)
	# capture history matrix
	write(t(chmat),con,ncolumns=nocc,append=TRUE)
	# first occasions seen 
	write(model_data$imat$first,con,ncolumns=n,append=TRUE)
	# frequency of capture history 
	if(!re)
	{
		write(model_data$imat$freq,con,ncolumns=n,append=TRUE)
	} else
	{
		if(any(model_data$imat$freq!=1))stop("\n cannot use random effects with frequency >1")
	}
	write(t(model_data$time.intervals),con,ncolumns=nocc-1,append=TRUE)
	# S design matrix
	write(ncol(model_data$S.dm),con,append=TRUE)
	write(t(model_data$S.dm),con,ncolumns=ncol(model_data$S.dm),append=TRUE)
	# p design matrix
    # zero out dm if any unobserved stratum; only done to remove unneeded columns
    if(uS>0)
	{
		model_data$p.dm[as.numeric(ddl$p$stratum)>=(length(strata.labels)-uS),]=0
		select=vector("logical",length=ncol(model_data$p.dm))
		for (i in 1:ncol(model_data$p.dm))
			select[i]=any(model_data$p.dm[,i]!=0)
		model_data$p.dm=model_data$p.dm[,select,drop=FALSE]
	}
    write(ncol(model_data$p.dm),con,append=TRUE)
	write(t(model_data$p.dm),con,ncolumns=ncol(model_data$p.dm),append=TRUE)
	# Psi design matrix
	# zero out subtracted stratum and remove any unneeded columns
    model_data$Psi.dm[!is.na(ddl$Psi$fix),]=0
	select=vector("logical",length=ncol(model_data$Psi.dm))
	for (i in 1:ncol(model_data$Psi.dm))
		select[i]=any(model_data$Psi.dm[,i]!=0)
	model_data$Psi.dm=model_data$Psi.dm[,select,drop=FALSE]
	write(ncol(model_data$Psi.dm),con,append=TRUE)
	write(t(model_data$Psi.dm),con,ncolumns=ncol(model_data$Psi.dm),append=TRUE)
	close(con)
	con=file(paste(tpl,".pin",sep=""),open="wt")
	write(par$S,con,ncolumns=length(par$S),append=FALSE)
	write(par$p,con,ncolumns=length(par$p),append=TRUE)
	write(par$Psi,con,ncolumns=length(par$Psi),append=FALSE)
	close(con)   
	if(hessian)
		xx=run_admb(tpl,extra.args=extra.args)
	else
		xx=run_admb(tpl,extra.args=paste(extra.args,"-nohess"))
	convergence=attr(xx,"status")
	if(is.null(convergence))convergence=0
	res=read_admb(tpl)
	beta=list(unscale.par(c(res$coeflist$phibeta,res$coeflist$pbeta,res$coeflist$psibeta),scale))
	parnames=names(unlist(beta))
	fixed.npar=length(unlist(beta))
	if(!is.null(res$hes))
	{
		beta.vcv=solvecov(res$hes)$inv
		rownames(res$hes)=parnames
		colnames(res$hes)=rownames(res$hes)
		if(all(diag(beta.vcv>0))) 
			res$cor=beta.vcv/outer(sqrt(diag(beta.vcv)),sqrt(diag(beta.vcv)))
	}  else
		beta.vcv=res$vcov
	rownames(beta.vcv)=parnames
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
	results$coeflist=res$coeflist
	results$npar=list(npar=res$npar,npar_sdrpt=res$npar_sdrpt,npar_total=res$npar_total)
	results$beta.vcv=beta.vcv
	res=results
#  Restore non-accumulated, non-scaled dm's etc
   res$model_data=model_data.save
#  Assign S3 class values and return
   class(res)=c("crm","admb","mle","mscjs")
   return(res)
}


