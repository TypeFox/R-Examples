#' Constructs and runs a set of MARK models from a dataframe of parameter
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
#' @param run if FALSE, only calls \code{\link{make.mark.model}} to test for
#' model errors but does not run the models
#' @param use.initial if TRUE, initial values are constructed for new models
#' using completed models that have already been run in the set
#' @param initial vector, mark model or marklist for defining initial values
#' @param ... arguments to be passed to \code{\link{mark}}. These must be
#' specified as argument=value pairs.
#' @return if(run) marklist - list of mark models that were run and a model.table of
#' results; if(!run) a list of models that were constructed but not run.
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{collect.models}}, \code{\link{mark}},
#' \code{\link{create.model.list}}
#' @keywords utility
mark.wrapper <-
function(model.list,silent=FALSE,run=TRUE,use.initial=FALSE,initial=NULL,...)
{
# -----------------------------------------------------------------------------------------------------------------------
# mark.wrapper  -  a wrapper for the mark function; it takes all the arguments and passes them onto mark
#
#  Value:
#
#  returns a list of mark models
#
# -----------------------------------------------------------------------------------------------------------------------
initiallist=NULL
if(class(initial)[1]=="marklist")
	if(nrow(initial$model.table)!=nrow(model.list))
		stop("marklist specified for initial argument does not contain same number of models")
	else
		initiallist=initial
model.names=rep(NA,nrow(model.list))
if (!run) list.of.models<-vector("list", nrow(model.list))
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
  if(run)
  {
     mymodel=try(mark(model.parameters=model.parameters,initial=initial,silent=silent,...),silent=silent)
  }
  else
  {
     mymodel=try(make.mark.model(parameters=model.parameters,initial=initial,...),silent=silent)
	 list.of.models[[i]]<-mymodel
  }
  if(class(mymodel)[1]!="try-error")
  {
    eval(parse(text=paste(model.name,"=mymodel")))
	model.names[i]=model.name
    if(!run)
    {
       message("\n Design matrix columns: ", dim(mymodel$design.matrix)[2],"\n")
       print(colnames(mymodel$design.matrix))
    }
  }
}
rm(mymodel)
#
# Return fitted MARK model object
#
rm(initial)
if(run)
	if(silent)
      return(suppressMessages(collect.models()))
    else
	  return(collect.models())
else
   return(list.of.models)
}


#' Load external model
#' 
#' Loads external model into workspace with name model
#' 
#' 
#' @param model name of MARK model object stored externally
#' @return None; side effect only to load object with name model (not name
#' given to model)
#' @author Jeff Laake
#' @keywords utility
load.model=function(model)
{ 
  if(is.character(model))
  {
    if(file.exists(model))
       load(model)
    else
       stop(paste("Cannot find file",model))
  }
return(model)
}
