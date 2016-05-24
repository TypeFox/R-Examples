#' Collect names of MARK model objects from list of R objects (internal
#' function)
#' 
#' Either names of all \code{mark} model objects (\code{type=NULL}) or names of
#' \code{mark} model objects of a specific type (\code{type}) are extracted
#' from a vector of R objects (\code{lx}) that was collected from the parent
#' environment (frame) of the function that calls \code{collect.model.names}.
#' Thus, it is two frames back (parent.frame(2)).
#' 
#' If \code{type=NULL} then the names of all objects of
#' \code{class(x)[1]="mark"} in \code{lx} are returned.  If \code{type} is
#' specified, then the names of all objects of \code{class(x)=c("mark",type)}
#' in \code{lx} are returned.
#' 
#' This function was written with the intention that it would be called from
#' other functions ( e.g., \code{\link{collect.models}},
#' \code{\link{run.models}}) but it will work if called directly (e.g.,
#' \code{collect.model.names( lx=ls())}). While this function returns a vector
#' of model names, \code{\link{collect.models}} returns a list of model
#' objects.  The latter can be used to easily create a list of models created
#' in a function to be used as a return value without listing all the names of
#' the functions.  It uses \code{collect.model.names} to perform that function.
#' 
#' @param lx vector of R object names from parent.frame(2)
#' @param type either NULL (for all types) or a character model type (eg "CJS")
#' @param warning if TRUE warning given when models of different types are
#' collected
#' @return model.list: a vector of \code{mark} model names
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{collect.models}}, \code{\link{run.models}},
#' \code{\link{model.table}}
#' @keywords utility
collect.model.names <-
function(lx, type=NULL, warning=TRUE)
# ----------------------------------------------------------------------------------------
#
# collect.model.names - this function collects mark models of a particular type (eg "CJS")
#                       or if type is NULL it collects all MARK model names from a list of 
#                       objects (lx) collected from the parent environment
#                       of the function that calls this one (ie parent.frame(2) two back from here).
#                  
# Arguments:        
#
#  lx             - listing of R objects
#  type           - either NULL(all models) or a specific model type (eg "CJS")
#  warning        - if TRUE, gives message that models not found
#
# Value:
#
#  model.list     - a list of MARK model names to use
#
# Functions used: setup.model
#
# ----------------------------------------------------------------------------------------
{
#
# Make sure type is correct if not NULL and assign to model
#
if(!is.null(type)) setup.model(type,0)
#
# Collect models from parent.frame of a particular "type" if any.
#
model.list=NULL
exclude=grep("\\*",lx)
if(length(exclude)>0)lx=lx[-exclude]
for (i in 1:length(lx)) 
{
   classval=class(eval(parse(text=lx[i]),envir=parent.frame(2)))
   if(classval[1]=="mark")
     if(is.null(type))
        model.list=c(model.list,lx[i])
     else
        if(classval[2]==type)
           model.list=c(model.list,lx[i])
}
#
# 09 Jan 06; Changed stops to cat statements and return NULL if none found
#
if(length(model.list)==0)
{
  if(warning)
  {
     if(is.null(type))
        message(paste("\nNo mark models found\n"))
     else
        message(paste("\nNo",type,"mark models found\n"))
  }
  return(NULL)
}
else
   return(model.list)
}
