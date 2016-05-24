#' Setup parameter structure specific to model (internal use)
#' 
#' Defines list of parameters used in the specified type of model
#' (\code{model}) and adds default values for each parameter to the list of
#' user specified values (eg formula, link etc).
#' 
#' @param model type of model ("CJS", "Burnham" etc)
#' @param parameters list of model parameter specifications
#' @param nocc number of occasions (value only specified if needed)
#' @param check if TRUE only the vector of parameter names is returned
#' \code{par.list}
#' @param number.of.groups number of groups defined for data
#' @export
#' @return The return value depends on the argument \code{check}. If it is TRUE
#' then the return value is a vector of the names of the parameters used in the
#' specified type of model. For example, if \code{model="CJS"} then the return
#' value is \code{c("Phi","p")}.  This is used by the function
#' \code{\link{valid.parameters}} to make sure that parameter specifications
#' are valid for the model (i.e., specifying recovery rate r for "CJS" would
#' give an error).  If the function is called with the default of
#' \code{check=FALSE}, the function returns a list of parameter specifications
#' which is a modification of the argument \code{parameters} which adds
#' parameters not specified and default values for all types of parameters that
#' were not specified. The list length and names of the list elements depends
#' on the type of model. Each element of the list is itself a list with varying
#' numbers of elements which depend on the type of parameter although some
#' elements are the same for all parameters.  Below the return value list is
#' shown generically with parameters named p1,...,pk.  \tabular{ll}{ \code{p1}
#' \tab List of specifications for parameter 1 \cr \code{p2} \tab List of
#' specifications for parameter 2 \cr . \tab \cr .  \tab \cr .  \tab \cr
#' \code{pk} \tab List of specifications for parameter k \cr }
#' 
#' The elements for each parameter list all include: \tabular{ll}{ \code{begin}
#' \tab 0 or 1; beginning time for the first \cr \tab parameter relative to
#' first occasion \cr \code{num} \tab 0 or -1; number of parameters relative to
#' \cr \tab number of occassions \cr \code{type} \tab type of PIM structure;
#' either "Triang" or "Square" \cr \code{formula} \tab formula for parameter
#' model (e.g., \code{~time}) \cr \code{link} \tab link function for parameter
#' (e.g., \code{"logit"}) \cr }
#' 
#' and may include: \tabular{ll}{ \code{share} \tab only valid for p in closed
#' capture models; \cr \tab if TRUE p and c models shared \cr \code{mix} \tab
#' only valid for closed capture heterogeneity \cr \tab models; if TRUE
#' mixtures are used \cr \code{rows} \tab only valid for closed capture
#' heterogeneity models \cr \code{fixed} \tab fixed values specified by user
#' and \cr \tab not used modified in this function \cr }
#' @author Jeff Laake
#' @seealso \code{\link{setup.model}},\code{\link{valid.parameters}}
#' @keywords utility
setup.parameters <-
		function(model,parameters=list(),nocc=NULL,check=FALSE,number.of.groups=1)
# ----------------------------------------------------------------------------------------
#  setup.parameters  - fills in value for begin and num for each parameter type depending
#                      on the type of c-r model. num defines number of parameters relative to
#                      number of occasions.  begin defines the first occasion which is relevant
#                      to the parameter
#
#  Arguments:
#    model      - type of model ("CJS", "JS" etc)
#    parameters - list of model parameter specifications
#    nocc       - number of occasions (value only specified if needed)
#    check      - default is FALSE; if TRUE it only returns list of parameter names
#
#  Value:
#    parameters - updated list of model parameter specifications with new fields added
#
#
# ----------------------------------------------------------------------------------------
{
# Read in parameter definitions
	fdir=system.file(package="marked")	
	fdir=file.path(fdir,"parameters.txt")	
	parameter_definitions=read.delim(fdir,header=TRUE,
			colClasses=c("character","character",rep("numeric",3),rep("character",3),
					rep("logical",3),"numeric","logical","logical","character","character","logical"))
#
#  Create valid parameter list depending on model.
#
	parameter_definitions=parameter_definitions[toupper(parameter_definitions$model)==model,]
	par.list=parameter_definitions$parname
#
#  If this is just a parameter check, return par.list
#
	if(check)return(par.list)
#
#  For each parameter create an empty list if none specified in input
#   
	pars=vector("list",length(par.list))
	names(pars)=par.list
	for (i in 1:length(par.list))
	{
		if(par.list[i]%in%names(parameters))
			pars[[i]]=parameters[[par.list[i]]]
	}
#
#  Next depending on model type, assign non-specified default values
#
	for(i in 1:length(par.list))
	{
		for(j in 3:ncol(parameter_definitions))		
			if(!is.na(parameter_definitions[i,j]) & parameter_definitions[i,j]!="" & !names(parameter_definitions)[j]%in%names(parameters[[par.list[i]]]))
				pars[[par.list[i]]][names(parameter_definitions)[j]]=list(parameter_definitions[i,j])
		if(pars[[par.list[i]]]$formula!=" " && is.character(pars[[par.list[i]]]$formula))
			pars[[par.list[i]]]$formula=as.formula(pars[[par.list[i]]]$formula)
		if(is.null(pars[[par.list[i]]]$num))pars[[par.list[i]]]$num=NA
		if(!is.na(pars[[par.list[i]]]$num)&&pars[[par.list[i]]]$num==1)pars[[par.list[i]]]$num=-(nocc-1)
		if(!is.null(pars[[par.list[i]]]$share) && pars[[par.list[i]]]$share && is.null(pars[[par.list[i]]]$pair)) pars[[par.list[i]]]$share=NULL
#
#       if include or mlogit have multiple values turn into a vector
#
		if(!is.null(pars[[par.list[i]]]$include))
		{
			xx=strsplit(pars[[par.list[i]]]$include,",")
			if(length(xx)==1)
				pars[[par.list[i]]]$include=xx[[1]]
		}
		if(!is.null(pars[[par.list[i]]]$mlogit))
		{
			xx=strsplit(pars[[par.list[i]]]$mlogit,",")
			if(length(xx)==1)
				pars[[par.list[i]]]$mlogit=xx[[1]]
		}
	}
	return(pars)
}
