#' Defines model specific parameters (internal use)
#' 
#' Compares \code{model}, the name of the type of model (eg "CJS") to the list
#' of acceptable models to determine if it is supported and then creates some
#' global fields specific to that type of model that are used to modify the
#' operation of the code.
#' 
#' In general, the structure of the different types of models (e.g.,
#' "CJS","Recovery",...etc) are very similar with some minor exceptions.  This
#' function is not intended to be called directly by the user but it is
#' documented to enable other models to be added.  This function is called by
#' other functions to validate and setup model specific parameters.  For
#' example, for live/dead models, the length of the capture history is twice
#' the number of capture occasions and the number of time intervals equals the
#' number of capture occasions because the final interval is included with dead
#' recoveries.  Whereas, for recapture models, the length of the capture
#' history is the number of capture occasions and the number of time intervals
#' is 1 less than the number of occasions.  This function validates that the
#' model is valid and sets up some parameters specific to the model that are
#' used in the code.
#' 
#' @param model name of model type (must be in vector \code{valid.models})
#' @param nocc length of capture history string
#' @param mixtures number of mixtures
#' @export
#' @aliases setup.model setupHMM
#' @return model.list - a list with following elements \item{etype}{encounter
#' type string for MARK input; typically same as model} \item{nocc}{number of
#' capture occasions} \item{num}{number of time intervals relative to number of
#' occasions (0 or -1)} \item{mixtures}{number of mixtures if any}
#' \item{derived}{logical; TRUE if model produces derived estimates}
#' @author Jeff Laake
#' @seealso \code{\link{setup.parameters}}, \code{\link{valid.parameters}}
#' @keywords utility
setup.model <-
		function(model,nocc,mixtures=1)
{
    # Read in parameter definitions
	fdir=system.file(package="marked")	
	fdir=file.path(fdir,"models.txt")	
	model_definitions=read.delim(fdir,header=TRUE,
			colClasses=c("character",rep("numeric",1),rep("logical",4)))
	model_def=model_definitions[model_definitions$model==model,]	
	if(nrow(model_def)==0)
		stop("Invalid type of model = ",model," Valid types are\n", paste(model_definitions$model,collapse="\n"))
    # model_def$nocc=nocc/model_def$divisor; not used at present
	model_def$nocc=nocc
	model_def=as.list(model_def)
	return(model_def)
}
setupHMM=function(model_def,model,strata.labels)
{
	if(toupper(model)=="HMMCJS")
	{
		model_def$hmm$ObsLevels=c(0,1)
		model_def$hmm$m=2
		model_def$hmm$fct_dmat=cjs_dmat
		model_def$hmm$fct_gamma=cjs_gamma
		model_def$hmm$fct_delta=cjs_delta
	} 
	if(toupper(model)%in%c("HMMCJS1TL"))
	{
		model_def$hmm$ObsLevels=c("+","-","0")
		model_def$hmm$fct_dmat=cjs1tl_dmat
		model_def$hmm$fct_gamma=cjs1tl_gamma
		model_def$hmm$fct_delta=cjs_delta
		model_def$hmm$strata.labels=c("1","0")
		model_def$hmm$obs_strata_map=1:2
		model_def$hmm$m=3
	}
	if(toupper(model)%in%c("HMMCJS2TL"))
	{
		model_def$hmm$ObsLevels=c("++","+-","-+","--","0")
		model_def$hmm$fct_dmat=cjs2tl_dmat
		model_def$hmm$fct_gamma=cjs2tl_gamma
		model_def$hmm$fct_delta=cjs_delta
		model_def$hmm$strata.labels=c("11","10","01","00")
		model_def$hmm$strata_data=data.frame(tag1=c(0,0,1,1),tag2=c(0,1,0,1))
		model_def$hmm$obs_strata_map=1:4
		model_def$hmm$m=5
	}
	if(toupper(model)%in%c("BAYESMSCJS","HMMMSCJS"))
	{
		model_def$hmm$ObsLevels=c(0,strata.labels)
		model_def$hmm$fct_dmat=ms_dmat
		model_def$hmm$fct_gamma=ms_gamma
		model_def$hmm$fct_delta=cjs_delta
		model_def$hmm$strata.labels=strata.labels
		model_def$hmm$m=length(strata.labels)+1
	}   
	if(toupper(model)=="HMMUMSCJS")
	{
		model_def$hmm$ObsLevels=c(0,strata.labels,"U")
		model_def$hmm$fct_dmat=ums_dmat
		model_def$hmm$fct_gamma=ms_gamma
		model_def$hmm$fct_delta=cjs_delta
		model_def$hmm$strata.labels=strata.labels
		model_def$hmm$m=length(strata.labels)+1
	}
	if(toupper(model)%in%c("HMMU2MSCJS","HMMU2IMSCJS"))
	{
		if(length(strata.labels)!=2 | !"states"%in%names(strata.labels))
			stop("structure of strata labels is incorrect; list of 2 character vectors with one named states")
		obs.states=c(strata.labels$states,"U")
		strata=strata.labels[[names(strata.labels)[names(strata.labels)!="states"]]]
        model_def$hmm$ObsLevels=c(0,apply(rev(expand.grid(list(obs.states,strata))),1,paste,collapse=""))
		model_def$hmm$fct_delta=cjs_delta		
		model_def$hmm$fct_dmat=ums2_dmat		
		model_def$hmm$fct_gamma=ms_gamma
		model_def$hmm$strata.list=strata.labels
		model_def$hmm$strata.labels=apply(rev(expand.grid(list(strata.labels$states,strata))),1,paste,collapse="")
		model_def$hmm$m=length(model_def$hmm$strata.labels)+1
	}
	if(toupper(model)=="HMMU2IMSCJS") model_def$hmm$fct_gamma=ms2_gamma
	if(toupper(model)=="MVMSCJS")
	{
		model_def$hmm$fct_dmat=mvms_dmat
		model_def$hmm$fct_gamma=mvms_gamma
		model_def$hmm$fct_delta=cjs_delta
		model_def$hmm$fct_sup=mvms_sup
		model_def$hmm$strata.list=set_mvms(strata.labels)
		model_def$hmm$strata.labels=apply(model_def$hmm$strata.list$df.states,1,paste,collapse="")
		model_def$hmm$m=nrow(model_def$hmm$strata.list$df.states)+1
		model_def$hmm$ObsLevels=c(0,apply(model_def$hmm$strata.list$df,1,paste,collapse=""))
	}
	if(toupper(model)=="MVMS")
	{
		model_def$hmm$fct_dmat=mvms_dmat
		model_def$hmm$fct_gamma=mvms_gamma
		#model_def$hmm$fct_delta=pi_mat
		model_def$hmm$strata.list=set_mvms(strata.labels)
		model_def$hmm$strata.labels=apply(model_def$hmm$strata.list$df.states,1,paste,collapse="")
		model_def$hmm$m=nrow(model_def$hmm$strata.list$df.states)+1
		model_def$hmm$ObsLevels=c(0,apply(model_def$hmm$strata.list$df,1,paste,collapse=""))
	}
	return(model_def)
}



