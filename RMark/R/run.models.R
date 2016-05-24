#' Runs a set of MARK models
#' 
#' Runs either a collection of models as defined in \code{model.list} or runs
#' all defined MARK object models in the frame of the calling function with no
#' \code{output} (\code{model.list=NULL}) or just those of a particular type
#' (e.g., \code{type="CJS"})
#' 
#' The model names in \code{model.list} must be in the frame of the function
#' that calls \code{run.models}. If \code{model.list=NULL} or the MARK models
#' are collected from the frame of the calling function (the parent). If
#' \code{type} is specified only the models of that type (e.g., "CJS") are run.
#' In each case the models are run and saved in the parent frame.
#' 
#' @param model.list either a vector of model names or NULL to run all MARK
#' models possibly of a particular \code{type}
#' @param type either a model type (eg "CJS", "Burnham" or "Barker") or NULL
#' for all types
#' @param save if TRUE, the R data directory is saved (i.e.,
#' \code{save.image()}) between analyses to enable interruption without losing
#' analyses that have already been run
#' @param ... any additional parameters to be passed to
#' \code{\link{run.mark.model}}
#' @return None; models are stored in parent frame.
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{collect.model.names}}, \code{\link{run.mark.model}}
#' @keywords utility
run.models <-
function(model.list=NULL,type=NULL, save=TRUE, ...)
# -----------------------------------------------------------------------------------------------------------------------
#
# run.models   - runs either a list of models as defined in model.list or runs all
#                defined models with no $output defined for all MARK objects (model.list=NULL)
#                of just those of a particular type (eg type="CJS")
#
# The model objects are run and saved in the calling environment. If save=TRUE, a save.image() is
# done between runs in case there is a problem.
#
# Arguments:
#
#   model.list  - either a vector of model names to run or NULL to run all MARK models
#   type        - either a model type (eg "CJS", "Burnham" or "Barker") or NULL for all types
#   save        - if TRUE a save.image() is done between analyses
#   ...         - additonal parameters for call to mark
#
# Value: 
#   None
#
# Functions used: collect.model.names, run.mark.model  
#
# -----------------------------------------------------------------------------------------------------------------------
{
#
# Get a list of objects from calling frame
#  
  lx=ls(envir=parent.frame())
#
# Collect appropriate MARK models from frame list unless a list has already been given
#
  model.list=collect.model.names(lx, type)
  if(is.null(model.list)) stop("No models need to be run")
  run=FALSE
#
#  For each model in the list, run the model unless it already contains output
#
  for(i in 1:length(model.list))
  {
     model=eval(parse(text=model.list[i]),envir=parent.frame())
     if(is.null(model$output))
     {
        run=TRUE
        eval(parse(text=paste(model.list[i],"=run.mark.model(",model.list[i],"...)")),envir=parent.frame())
        if(save)save.image()
     }
  }
  if(!run)message("All defined models have been run\n")
invisible()
}
