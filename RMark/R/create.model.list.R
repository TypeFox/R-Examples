#' Creates a dataframe of all combinations of parameter specifications
#' 
#' Creates a dataframe of all combinations of parameter specifications for each
#' parameter in a particular type of MARK \code{model}. It is used together
#' with \code{\link{mark.wrapper}} to run a series of models from sets of
#' parameter specifications.
#' 
#' This function scans the frame of the calling enviroment and collects all
#' list objects that contain a formula and have names that match
#' \code{parameter.} where parameter is the name of a type of parameter in the
#' \code{model} type.  For example, it looks for \code{Phi.} and \code{p.} for
#' \code{model="CJS"}. Any number of characters can follow the period.  Each of
#' the named objects should specify a list that matches the structure of a
#' parameter specification as described in \code{\link{make.mark.model}}. It
#' only collects list objects that contain an element named \code{formula},
#' thus it will not collect one like \code{Phi.fixed=list(fixed=1)}.  If you
#' want to do something like that, specify it as
#' \code{Phi.fixed=list(formula=~1,fixed=1)}. It is safest to use this inside a
#' function that defines all of the parameter specifications as shown in the
#' example below.  The primary use for this function is to create a dataframe
#' which is passed to \code{\link{mark.wrapper}} to construct and run each of
#' the models.  It was written as a separate function to provide flexibility to
#' add/delete/modify the list prior to passing to \code{\link{mark.wrapper}}.
#' For example, only certain combinations may make sense for some parameter
#' specifications.  Thus you could define a set to create all the combinations
#' and then delete the ones from the dataframe that do not make sense. you
#' want, add others and re-run the function and merge the resulting dataframes.
#' If there are no specifications found for a particular model parameter, it is
#' not included in the list and when it is passed to
#' \code{\link{make.mark.model}}, the default specification will be used for
#' that parameter.
#' 
#' @param model character string identifying the type of model (e.g., "CJS")
#' @return dataframe of all combinations of parameter specifications for a
#' model.  Each field (column) is the name of a type of parameter (e.g., p and
#' Phi for CJS).  The values are character strings identifying particular
#' parameter specifications.
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{mark.wrapper}}
#' @keywords utility
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' #
#' # Compare this to the run.dipper shown under ?dipper
#' # It is only necessary to create each parameter specification and 
#' # create.model.list and mark.wrapper will create and run models for 
#' # each combination. Notice that the naming of the parameter 
#' # specifications has been changed to accommodate format for 
#' # create.model.list. Only a subset of the parameter specifications 
#' # are used here in comparison to other run.dipper
#' #
#' data(dipper)
#' run.dipper=function()
#' {
#'      #
#'      # Process data
#'      #
#'      dipper.processed=process.data(dipper,groups=("sex"))
#'      #
#'      # Create default design data
#'      #
#'      dipper.ddl=make.design.data(dipper.processed)
#'      #
#'      # Add Flood covariates for Phi and p that have different values
#'      #
#'      dipper.ddl$Phi$Flood=0
#'      dipper.ddl$Phi$Flood[dipper.ddl$Phi$time==2 |
#'                   dipper.ddl$Phi$time==3]=1
#'      dipper.ddl$p$Flood=0
#'      dipper.ddl$p$Flood[dipper.ddl$p$time==3]=1
#'      #
#'      #  Define range of models for Phi
#'      #
#'      Phi.dot=list(formula=~1)
#'      Phi.time=list(formula=~time)
#'      Phi.sex=list(formula=~sex)
#'      Phi.Flood=list(formula=~Flood)
#'      #
#'      #  Define range of models for p
#'      #
#'      p.dot=list(formula=~1)
#'      p.time=list(formula=~time)
#'      p.Flood=list(formula=~Flood)
#'      #
#'      # Run all pairings of models
#'      #
#'      dipper.model.list=create.model.list("CJS")
#'      dipper.results=mark.wrapper(dipper.model.list,
#'               data=dipper.processed,ddl=dipper.ddl)
#'      #
#'      # Return model table and list of models
#'      #
#'      return(dipper.results)
#' }
#' dipper.results=run.dipper()
#' }
create.model.list<-function(model)
{
	parameters=setup.parameters(model,check=TRUE)
	model.list=list()
	for(n in parameters)
	{
		vec=ls(pattern=paste("^",n,"\\.",sep=""),envir=parent.frame())
		if(length(vec)>0)
		{
			for (i in 1:length(vec))
			{
				if(eval(parse(text=paste("is.list(",vec[i],")",sep="")),envir=parent.frame()))
				{
					if(eval(parse(text=paste("!is.null(",vec[i],"$formula)",sep="")),envir=parent.frame()) |
							(eval(parse(text=paste("is.list(",vec[i],"[[1]])",sep="")),envir=parent.frame())&&
							eval(parse(text=paste("!is.null(",vec[i],"[[1]]$formula)",sep="")),envir=parent.frame())))
						model.list[[n]]=c(model.list[[n]],vec[i])
				}
			}
		} else
			message("Using default formula for ",n,"\n")
	}
	if(length(model.list)==0)
		stop("\nNo model specifications found. Use case sensitive parameter.description notation (e.g., Phi.time)\n")
	if(length(model.list)>1)
	{
		model.list=expand.grid(model.list)
		for (j in 1:dim(model.list)[2])
			model.list[,j]=as.character(model.list[,j])
	}
	else
		model.list=as.data.frame(model.list)
	return(model.list)
}

