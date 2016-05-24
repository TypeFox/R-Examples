#' Constructs and runs in parallel a set of MARK models from a dataframe of parameter
#' specifications
#' 
#' This is a convenience function that uses a dataframe of parameter
#' specifications created by \code{\link{create.model.list}} and it constructs
#' and runs each model and names the models by concatenating each of the
#' parameter specification names separated by a period.  The results are
#' returned as a marklist with a model.table constructed by
#' \code{\link{collect.models}}.
#' 
#' The model names in \code{model.list} must be in the frame of the function
#' that calls \code{run.models}. If \code{model.list=NULL} or the MARK models
#' are collected from the frame of the calling function (the parent). If
#' \code{type} is specified only the models of that type (e.g., "CJS") are run.
#' In each case the models are run and saved in the parent frame. To fully
#' understand, how this function works and its limitations, see
#' \code{\link{create.model.list}}.
#' 
#' If \code{use.initial=TRUE}, prior to running a model it looks for the first
#' model that has already been run (if any) for each parameter formula and
#' constructs an \code{initial} vector from that previous run. For example, if
#' you provided 5 models for p and 3 for Phi in a CJS model, as soon as the
#' first model for p is run, in the subsequent 2 models with different Phi
#' models, the initial values for p are assigned based on the run with the
#' first Phi model.  At the outset this seemed like a good idea to speed up
#' execution times, but from the one set of examples I ran where several
#' parameters were at boundaries, the results were discouraging because the
#' models converged to a sub-optimal likelihood value than the runs using the
#' default initial values.  I've left sthis option in but set its default value
#' to FALSE.
#' 
#' A possibly more useful argument is the argument \code{initial}.  Previously,
#' you could use \code{initial=model} as part of the ... arguments and it would
#' use the estimates from that model to assign initial values for any model in
#' the set. Now I've defined \code{initial} as a specific argument and it can
#' be used as above or you can also use it to specify a \code{marklist} of
#' previously run models. When you do that, the code will lookup each new model
#' to be run in the set of models specified by \code{initial} and if it finds
#' one with the matching name then it will use the estimates for any matching
#' parameters as initial values in the same way as \code{initial=model} does.
#' The model name is based on concatenating the names of each of the parameter
#' specification objects.  To make this useful, you'll want to adapt to an
#' approach that I've started to use of naming the objects something like
#' p.1,p.2 etc rather than naming them something like p.dot, p.time as done in
#' many of the examples.  I've found that using numeric approach is much less
#' typing and cumbersome rather than trying to reflect the formula in the name.
#' By default, the formula is shown in the model selection results table, so it
#' was a bit redundant.  Now where I see this being the most benefit.
#' Individual covariate models tend to run rather slowly. So one approach is to
#' run the sequence of models (eg results stored in initial_marklist),
#' including the set of formulas with all of the variables other than
#' individual covariates.  Then run another set with the same numbering scheme,
#' but adding the individual covariates to the formula and using
#' \code{initial=initial_marklist} That will work if each parameter
#' specification has the same name (eg., p.1=list(formula=~time) and then
#' p.1=list(formula=~time+an_indiv_covariate)).  All of the initial values will
#' be assigned for the previous run except for any added parameters (eg.
#' an_indiv_covariate) which will start with a 0 initial value.
#' 
#' @param model.list a dataframe of parameter specification names in the
#' calling frame
#' @param silent if TRUE, errors that are encountered are suppressed
#' @param use.initial if TRUE, initial values are constructed for new models
#' using completed models that have already been run in the set
#' @param initial vector, mark model or marklist for defining initial values
#' @param parallel if TRUE, runs models in parallel on multiple cpus
#' @param cpus number of cpus to use in parallel
#' @param threads number of cpus to use with mark.exe if positive or number of cpus to remain idle if negative
#' @param ... arguments to be passed to \code{\link{mark}}. These must be
#' specified as argument=value pairs.
#' @return marklist - list of mark models that were run and a model.table of
#' results
#' @author Eldar Rakhimberdiev
#' @export 
#' @import snowfall
#' @seealso \code{\link{collect.models}}, \code{\link{mark}},
#' \code{\link{create.model.list}}
#' @keywords utility
#' @examples 
#' \donttest{
#' # example not run to reduce time required for checking
#' do.MSOccupancy=function()
#' {
#' #  Get the data
#' 	data(NicholsMSOccupancy)
#' # Define the models; default of Psi1=~1 and Psi2=~1 is assumed
#' # p varies by time but p1t=p2t
#' 	p1.p2equal.by.time=list(formula=~time,share=TRUE)
#' # time-varying p1t and p2t
#' 	p1.p2.different.time=list(p1=list(formula=~time,share=FALSE),p2=list(formula=~time))
#' #  delta2 model with one rate for times 1-2 and another for times 3-5;
#' # delta2 defined below
#' 	Delta.delta2=list(formula=~delta2)
#' 	Delta.dot=list(formula=~1)  # constant delta
#' 	Delta.time=list(formula=~time) # time-varying delta
#' # Process the data for the MSOccupancy model
#' 	NicholsMS.proc=process.data(NicholsMSOccupancy,model="MSOccupancy")
#' # Create the default design data
#' 	NicholsMS.ddl=make.design.data(NicholsMS.proc)
#' # Add a field for the Delta design data called delta2.  It is a factor variable
#' # with 2 levels: times 1-2, and times 3-5.
#' 	NicholsMS.ddl=add.design.data(NicholsMS.proc,NicholsMS.ddl,"Delta",
#' 			type="time",bins=c(0,2,5),name="delta2")
#' # Create a list using the 4 p modls and 3 delta models (12 models total)
#' 	cml=create.model.list("MSOccupancy")
#' # Fit each model in the list and return the results
#' 	return(mark.wrapper.parallel(cml,data=NicholsMS.proc,ddl=NicholsMS.ddl,
#'     cpus=2,parallel=TRUE))
#' }
#' xx=do.MSOccupancy()
#' }
mark.wrapper.parallel<-
		function(model.list,silent=FALSE,use.initial=FALSE,initial=NULL, parallel=TRUE, cpus=2, threads=1, ...)
{
	
# -----------------------------------------------------------------------------------------------------------------------
# mark.wrapper  -  a wrapper for the mark function; it takes all the arguments and passes them onto mark
#
#  Value:
#
#  returns a list of mark models
#
# -----------------------------------------------------------------------------------------------------------------------
	if(R.Version()$os!="mingw32")
	{
		cat("\nWindows only function. Unable to get this function to run on non-Windows machine\n")
		return(NULL)
	}
	args<-match.call()
	t.1<-which(as.character(args[-1]) %in% ls(envir=parent.frame()))
	t.2<-mget(as.character(args[-1][t.1]), envir=parent.frame())
	names(t.2)<-names(args[-1][t.1])
	if (length(t.2)==length(args[-1])) {list.args<-t.2}
	else {
		list.args<-c(t.2, as.list(args[-1][!names(args[-1]) %in% names(t.2)]))
	}
	if (!any(names(list.args)=="parallel")) list.args$parallel=TRUE
	if (!any(names(list.args)=="wd")) list.args$wd=getwd()
	if (!any(names(list.args)=="prefix")) list.args$prefix="mark"
	existing.files<-list.files(path=list.args$wd, pattern=paste(list.args$prefix, "[[:digit:]]+\\.", sep=""))
	Max.number<-0
	if (!length(existing.files)==0) {
		Max.number<-sum(0, max(as.numeric(sapply(strsplit(substr(existing.files, 5, 100), "\\."), "[[", i=1)), na.rm=TRUE), na.rm=TRUE)
		if(list.args$prefix!="mark") warning("you already have files with prefix ", list.args$prefix, " in directory ", list.args$wd, ".\n first file will be named as ", list.args$prefix, formatC((Max.number+1), width=3,digits=0,format="f", flag="0"), "\n")
	}
	if (!any(names(list.args)=="silent")) list.args$silent=FALSE
	if (!any(names(list.args)=="run")) list.args$run=TRUE
	if (!any(names(list.args)=="use.initial")) list.args[-which(names(list.args)=="use.initial")]
	if (!any(names(list.args)=="invisible")) list.args$invisible=FALSE
	if (any(names(list.args)=="model.list")) list.args<-list.args[-which(names(list.args)=="model.list")]
	if (any(names(list.args)=="cpus")) list.args<-list.args[-which(names(list.args)=="cpus")]
	if (!any(names(list.args)=="threads")) list.args$threads<-1
#	if (any(names(list.args)=="threads")) list.args<-list.args[-which(names(list.args)=="threads")]
#	require("parallel")
	if (threads*cpus>parallel::detectCores() | (threads<0 & cpus!=1)) 
		stop("you've tried to use more cores than you have, try to combine threads and cpus to make ", 
				 parallel::detectCores(), " or less as a product\n")
	initiallist=NULL
	if(class(initial)[1]=="marklist")
		if(nrow(initial$model.table)!=nrow(model.list))
			stop("marklist specified for initial argument does not contain same number of models")
		else
			initiallist=initial
	model.names=rep(NA,nrow(model.list))
# for apply like fuctions we need list of models
	list.of.model.par<-vector("list", nrow(model.list))
	for (i in 1:nrow(model.list))
	{
		model.parameters=list()
		for(j in 1:ncol(model.list))
		{
			if(!is.list(eval(parse(text=model.list[i,j]),envir=parent.frame())[[1]]))
				model.parameters[[names(model.list)[j]]]=eval(parse(text=(as.character(model.list[i,j]))),envir=parent.frame())
		}
		for(j in 1:ncol(model.list))
		{
			if(is.list(eval(parse(text=model.list[i,j]),envir=parent.frame())[[1]]))
				model.parameters=c(model.parameters,eval(parse(text=(as.character(model.list[i,j]))),envir=parent.frame()))
		}
		model.name=paste(model.list[i,],collapse=".")
		if(!silent)message("\n",model.name,"\n")
		if(use.initial)
		{
			initial=NULL
			for(j in 1:ncol(model.list))
			{
				mindex=match(model.list[i,j],model.list[,j])
				if(!is.na(mindex)&& !is.na(model.names[mindex]))
				{
					estimates=eval(parse(text=paste(model.names[mindex],"$results$beta",sep="")))
					estimates=estimates[grep(paste(colnames(model.list)[j],":",sep=""),rownames(estimates)),]
					beta=estimates$estimate
					names(beta)=rownames(estimates)
					initial=c(initial,beta)
				}
			}
		}
		else
		if(!is.null(initiallist)) 
			if(model.name%in%names(initiallist))initial=initiallist[[model.name]]
			else
				initial=NULL
		if (any(names(list.args)=="initial")) list.args<-list.args[-which(names(list.args)=="initial")]
		list.of.model.par[[i]]<-c(list(model.parameters=model.parameters, initial=initial, model.name=NULL, filename=paste(list.args$prefix, formatC((Max.number+i), width=3,digits=0,format="f", flag="0"), sep="")), list.args) # here I add all of the args from wrapper
	}
	if (parallel) {
		sfInit(parallel=TRUE, cpus=cpus)
		suppressWarnings(sfLibrary("RMark",character.only=TRUE))
		list.of.models<-sfClusterApplyLB(list.of.model.par,  make.run.mark.model.apply.int)
		sfStop()
	}
	else list.of.models<-lapply(list.of.model.par,  make.run.mark.model.apply.int)
	rm(initial)
	#### here I am excluding all possible errors
	list.of.models<-list.of.models[!lapply(list.of.models, "[[", 1)=="error"]
	if(length(list.of.models)==0) stop("all models failed to run\n")
	modellist<-list.of.models[[1]]
	for(i in 2:(length(list.of.models))){
		tempobj=list.of.models[[i]]
		modellist<-merge.mark(modellist, tempobj)
		names(modellist)[i]<-paste("model",i, sep="." )
	}
	names(modellist)[1]<-"model.1"
	return(modellist)
}


make.run.mark.model.apply.int<-function(x) {
	old.wd<-getwd()
# here
	#require(RMark)
	run=as.logical(as.character(x[["run"]]))
	x<-x[-which(names(x)=="run")]
	if (any(names(x)=="silent")) {
		silent=as.logical(as.character(x[["silent"]]))
		x<-x[-which(names(x)=="silent")]
	}
	wd=x[["wd"]]
	x<-x[-which(names(x)=="wd")]
	setwd(as.character(wd))
	parallel=as.logical(as.character(x[["parallel"]]))
	x<-x[-which(names(x)=="parallel")]	
# removed code for !run which doesn't work and not needed.
#	if(run) {
		mymodel<-try(do.call(mark, x), silent=silent)
		if(class(mymodel)[1]=="try-error") mymodel<-"error"
		if (length(mymodel)>1)
			if(is.null(mymodel$results))
				mymodel<-"error"
#	}
#	else {
#		x<-x[-which(names(x)=="threads")]	
#		names(x)[names(x)=="model.parameters"]<-"parameters"
#		mymodel<-try(do.call(make.mark.model, x),silent=silent) }
	setwd(old.wd)
	return(mymodel)
}
