#' Automation of model runs
#' 
#' Some functions that help automate running a set of crm models based on parameter
#' specifications.
#' 
#' create.model.list creates all combinations of model specifications for the specified
#' set of parameters.  In the calling environment it looks for objects named parameter.xxxxxx where xxxxxx can
#' be anything. It creates a matrix with a column for each parameter and as many rows 
#' needed to create all combinations. This can be used as input to crm.wrapper.
#' 
#' crm.wrapper runs a sequence of crm models by constructing the call with the arguments
#' and the parameter specifications.  The parameter specifications can either be in the
#' local environment or in the environment of the named function models. The advantage of the
#' latter is that it is self-contained such that sets of parameter specifications can
#' be selected without possibility of being over-written or accidentally changed whereas 
#' with the former the set must be identified via a script and any in the environment will
#' be used which requires removing/recreating the set to be used.
#' 
#' @aliases  crm.wrapper create.model.list model.table load.model rerun_crm crmlist_fromfiles
#' @usage  crm.wrapper(model.list,data,ddl=NULL,models=NULL,base="",
#'               external=TRUE,run=TRUE,env=NULL,...)
#' 
#'         create.model.list(parameters)
#' 
#'         model.table(model.list)
#' 
#'         load.model(x)
#'  
#'         crmlist_fromfiles(filenames=NULL,external=TRUE)
#' 
#'         rerun_crm(data,ddl,model.list,method=NULL,modelnums=NULL,initial=NULL,...)
#' 
#' @param data Either the raw data which is a dataframe with at least one
#' column named ch (a character field containing the capture history) or a
#' processed dataframe. For rerun_crm this should be the processed dataframe
#' @param ddl Design data list which contains a list element for each parameter
#' type; if NULL it is created; For rerun_crm, must be the same ddl as used with original run can cannot be NULL
#' @param model.list matrix of model names contained in the environment of models function; each row is a model and each column is for a parameter and the value is formula name
#' @param models a function with a defined environment with model specifications as variables; values of model.list are some or all of those variables
#' @param base base value for model names
#' @param external if TRUE, model results are stored externally; otherwise they are stored in crmlist
#' @param run if TRUE, fit models; otherwise just create dml to test if model data are correct for formula
#' @param env environment to find model specifications if not parent.frame
#' @param ... aditional arguments passed to crm; for rerun_crm can be used to set hessian=TRUE for specific models after they have been run
#' @param parameters character vector of parameter names
#' @param x filename of externally stored model
#' @param method vector of methods to use for optimization if different that previous run in rerun_crm
#' @param modelnums model numbers to be re-run instead of those that did not covnerge
#' @param initial either a fitted crm model or the model number in model.list to use for starting values
#' @param filenames for non-Windows machine, vector of filenames for external files must be specifed in crmlist_fromfiles including .rda extension
#' @return create.model.list returns a matrix for crm.wrapper; crm.wrapper runs and stores models externally and retrurns a list of model results
#' and a model selection table; load.model returns model object that is stored externally
#' @author Jeff Laake
#' @export create.model.list
#' @export crm.wrapper
#' @export model.table
#' @export load.model
#' @export crmlist_fromfiles
#' @export rerun_crm
#' @import parallel
#' @seealso \code{\link{crm}}
#' @keywords models
crm.wrapper <- function(model.list,data,ddl=NULL,models=NULL,base="",external=TRUE,run=TRUE,env=NULL,...)
{
	if(is.null(env))env=parent.frame()
	results=vector("list",length=nrow(model.list)+1)
	results.names=NULL
	model.table=NULL
	for (i in 1:nrow(model.list))
	{
		model.parameters=list()		
		if(is.null(models))
		{
			for(j in 1:ncol(model.list))
			{
				if(!is.list(eval(parse(text=model.list[i,j]),envir=env)[[1]]))
					model.parameters[[names(model.list)[j]]]=eval(parse(text=(as.character(model.list[i,j]))),envir=env)
			}
			for(j in 1:ncol(model.list))
			{
				if(is.list(eval(parse(text=model.list[i,j]),envir=env)[[1]]))
					model.parameters=c(model.parameters,eval(parse(text=(as.character(model.list[i,j]))),envir=env))
			}	
		} else
		{
			model.parameters=models(model.list[i,])$model.parameters
		}
		model.name=paste(model.list[i,],collapse=".")
		cat(model.name,"\n")
		if(file.exists(paste(model.name,".rda",sep=""))&is.null(list(...)$initial))
		{
			load(paste(model.name,".rda",sep=""))
			initial=eval(parse(text=model.name))
			mymodel=crm(data=data,ddl=ddl,model.parameters=model.parameters,run=run,initial=initial,...)
		} else
		{
			mymodel=crm(data=data,ddl=ddl,model.parameters=model.parameters,run=run,...)
		}
		if(external)
		{
			assign(as.character(as.name(model.name)),mymodel)
			eval(parse(text=paste("save(",model.name,', file="',base,model.name,'.rda")',sep="")))
			results[[i]]=paste(model.name,".rda",sep="")
		} else
			results[[i]]=mymodel
		results.names=c(results.names,paste(model.list[i,],collapse="."))
		formulae=sapply(model.parameters,function(x){return(paste(x$formula,collapse=""))})
		formulae=paste(paste(names(formulae),"(",formulae,")",sep=""),collapse="")
		if(run)
		{
			df=data.frame(model=formulae,npar=sum(sapply(mymodel$results$beta,length)),AIC=mymodel$results$AIC,neg2lnl=mymodel$results$neg2lnl,convergence=mymodel$results$convergence)
			model.table=rbind(model.table,df)
		}
	}
	names(results)=results.names
	if(run)
	{
		model.table$DeltaAIC=model.table$AIC-min(model.table$AIC)
		model.table$weight=exp(-.5*model.table$DeltaAIC)
		model.table$weight=model.table$weight/sum(model.table$weight)
		model.table=model.table[order(model.table$DeltaAIC),c("model","npar","AIC","DeltaAIC","weight","neg2lnl","convergence")]
		results[[length(results)]]=model.table
		names(results)=c(results.names,"model.table")
		class(results)="crmlist"
	}
    return(results)
}
model.table=function(model.list=NULL)
{
    model.table=NULL
	for(i in 1:(length(model.list)-1))
	{
		if(!is.list(model.list[[i]]))
		{
			load(model.list[[i]])
			if(length(grep("\\\\",model.list[[i]]))>0 | length(grep("/",model.list[[i]]))>0)
			    model.list[[i]]=basename(model.list[[i]])	
			eval(parse(text=paste("mymodel=",unlist(strsplit(model.list[[i]],".rda")))))
			eval(parse(text = paste("rm(", unlist(strsplit(model.list[[i]], 
													".rda")),")")))
			gc()
		}else
			mymodel=model.list[[i]]
		formulae=sapply(mymodel$model.parameters,function(x){return(paste(x$formula,collapse=""))})
		formulae=paste(paste(names(formulae),"(",formulae,")",sep=""),collapse="")
		df=data.frame(model=formulae,npar=sum(sapply(mymodel$results$beta,length)),AIC=mymodel$results$AIC,neg2lnl=mymodel$results$neg2lnl,convergence=mymodel$results$convergence)
		model.table=rbind(model.table,df)
	}
	model.table$DeltaAIC=model.table$AIC-min(model.table$AIC)
	model.table$weight=exp(-.5*model.table$DeltaAIC)
	model.table$weight=model.table$weight/sum(model.table$weight)
	model.table=model.table[order(model.table$DeltaAIC),c("model","npar","AIC","DeltaAIC","weight","neg2lnl","convergence")]
	return(model.table)
}				
create.model.list<-function(parameters)
{
	model.list=list()
	for(n in parameters)
	{
		vec=ls(pattern=paste("^",n,"\\.",sep=""),envir=parent.frame())
		if(length(vec)>0)
			for (i in 1:length(vec))
			{
				if(eval(parse(text=paste("is.list(",vec[i],")",sep="")),envir=parent.frame()))
				{
					if(eval(parse(text=paste("!is.null(",vec[i],"$formula)",sep="")),envir=parent.frame()) |
							eval(parse(text=paste("!is.null(",vec[i],"[[1]]$formula)",sep="")),envir=parent.frame()))
						model.list[[n]]=c(model.list[[n]],vec[i])
				}
			}
		
	}
	if(length(model.list)==0)
		stop("\nNo model specifications found. Use parameter.description notation (e.g., Phi.time)\n")
	if(length(model.list)>1)
	{
		model.list=expand.grid(model.list)
		for (j in 1:dim(model.list)[2])
			model.list[,j]=as.character(model.list[,j])
	}
	else
		model.list=as.data.frame(model.list)
	model.list=model.list[ do.call(order, model.list) ,]
	return(model.list)
}
load.model=function(x)
{
	model=NULL
	if(!is.character(x)) stop("argument should be a character filename")
	if(!file.exists(x))stop(x,"does not exist in the working directory")
	load(x)
	model.name=strsplit(x,"\\.rda")[[1]][1]
	eval(parse(text=paste("assign(as.character(as.name('model')),",model.name,")")))
	return(model)
}
crmlist_fromfiles=function(filenames=NULL,external=TRUE)
{
	if(R.Version()$os!="mingw32")Filters=NULL
	if(is.null(filenames))
	{
		if(R.Version()$os=="mingw32")
			filenames=utils::choose.files(filters=Filters["RData",])
		else
			stop("filenames must be specified for non-Windows machine")
	}
	modelnames=unlist(strsplit(basename(filenames),"\\.rda"))
	mtable=model.table(c(filenames,""))
	model.list=vector("list",length=length(modelnames))
	names(model.list)=modelnames
	i=0
	for(f in filenames)
	{
		i=i+1
		m=modelnames[i]
		if(!external)
		{
  		    load(f)
			eval(parse(text=paste("model.list[['",m,"']]=",m,sep="")))
			eval(parse(text=paste("rm(",m,")",sep="")))
			gc()
		}
	    else
		   model.list[[m]]=m
	}
	model.list$model.table=mtable
	class(model.list)="crmlist"
	return(model.list)
}		
rerun_crm=function(data,ddl,model.list,method=NULL,modelnums=NULL,initial=NULL,...)
{
	# rerun through each model and if in modelnums (not NULL) or didn't 
	# converge (modelnums is NULL) rerun.
	for(i in 1:(length(model.list)-1))
	{
		if((is.null(modelnums) & model.list$model.table$convergence[i]!=0) | i %in% modelnums)
		{
			# if not a list then read in model
			model=model.list[[i]]
			external=FALSE
			if(!is.list(model)) 
			{
				model.name=unlist(strsplit(model.list[[i]],".rda"))
				external=TRUE
				load(model.list[[i]])
				eval(parse(text=paste("model=",model.name)))
			}
			# set default values
			if(is.null(model$results$mat))
				save.matrices=FALSE
			else
				save.matrices=TRUE
			if(is.null(method))
				method=attr(model$results$optim.details,"details")$method
			# set initial values to use; either current model or what is specified in initial
			# if external it is read in
			if(is.null(initial))
				initial=model
			else
			{
			   if(!is.list(initial)) 
			   {
				   if(is.numeric(initial))
				   {
					   if(length(initial)==1 && initial %in% 1:(length(model.list)-1))
						   initial=model.list[[initial]]
					   else
						   stop("Invalid initial value")
				   } else
				   {
					   imodel.name=unlist(strsplit(initial,".rda"))
					   load(initial)
					   eval(parse(text=paste("initial=",imodel.name)))
				   }
			   }
		    }
			# run model and save in mymodel
			mymodel=crm(data,ddl,model.parameters=model$model.parameters,design.parameters=model$design.parameters,initial=initial,
					save.matrices=save.matrices,method=method,...)
			# handle storing in list or externally saving
			if(external)
			{
				assign(as.character(as.name(model.name)),mymodel)
				eval(parse(text=paste("save(",model.name,', file="',model.name,'.rda")',sep="")))
				model.list[[i]]=paste(model.name,".rda",sep="")
			} else
				model.list[[i]]=mymodel
		}
	}
	model.list$model.table=model.table(model.list)
	return(model.list)
}
##' crm.wrapper.parallel can use multiple cpus to run multiple models simultaneously on nc cpus 
#crm.wrapper.parallel=function(nc,model.list,data,ddl=NULL,models=NULL,base="",external=TRUE,run=TRUE,...)
##' @param nc number of cpus to use for parallel model runs. Will not let you to use more than are available.
#{
#	split.model.list=split(model.list,1:nrow(model.list))
#	if(nc>detectCores())stop(paste("Too many cpus.  nc > ",detectCores(), " cpus available"))
#	cl=makeCluster(nc)
#	if(is.null(dim(cml)))
#		clusterExport(cl,unique(cml))
#	else
#	    clusterExport(cl,as.vector(apply(cml,2,unique)))
#	res<-parLapply(cl,split.model.list,crm.wrapper,data=data,ddl=ddl,models=models,base=base,external=external,run=run,...)
#	if(!external)
#	{
#		res=lapply(res,function(x) x[[1]])
#		res$model.table=model.table(res)
#	} else
#	{
#		res=paste(apply(model.list,1,paste,collapse="."),".rda",sep="")
#		names(res)=apply(model.list,1,paste,collapse=".")
#		res=as.list(res)
#		res$model.table=model.table(res)
#	}	
#	class(res)="crmlist"
#	stopCluster(cl)
#	return(res)
#}
##' @export crm.wrapper.parallel



