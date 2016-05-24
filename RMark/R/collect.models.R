#' Collect MARK models into a list and optionally construct a table of model
#' results
#' 
#' Collects \code{mark} models contained in \code{lx} of specified \code{type}
#' (if any) and returns models in a list with a table of model results if
#' \code{table=TRUE}.
#' 
#' If \code{lx} is NULL a vector of object names in the parent frame is
#' constructed for \code{lx}.  Within \code{lx} all \code{mark} model objects
#' (i.e., \code{class(x)[1]=="mark"}) are returned if \code{type} is NULL.  If
#' \code{type} is specified and is valid, then the names of all \code{mark}
#' model objects of the specified \code{type} (i.e., \code{class(x) =
#' c("mark",type)}) in \code{lx} are returned.  If \code{table=TRUE} a table of
#' model selection results is also included in the returned list.
#' 
#' This function was written to be able to easily collect a series of mark
#' models in a list without specifying the names of each model object.  This is
#' useful in constructing a return value of a function that runs a series of
#' models for a particular analysis. For an example see \code{\link{dipper}}.
#' 
#' @param lx if NULL, constructs vector of object names (\code{ls()}) from
#' frame of calling function otherwise it uses specifed names in \code{lx}
#' @param type either NULL (for all types) or a character model type (e.g.
#' \code{type="CJS"})
#' @param table if TRUE, a table of model results is also included in the
#' returned list
#' @param adjust if TRUE, adjusts number of parameters (npar) to number of
#' columns in design matrix which modifies AIC
#' @param external if TRUE the mark objects are saved externally rather than in
#' the resulting marklist; the filename for each object is kept in its place
#' @return model.list: a list of \code{mark} models and optionally a table of
#' model results.
#' @note This function and others that use it or use
#' \code{\link{collect.model.names}} to collect a series of models or assign a
#' value to a series of models (e.g., \code{\link{adjust.chat}}) should be used
#' with a degree of caution. It is important to understand the scope of the
#' collection.  If the call to this function is made at the R prompt, then it
#' will collect all models (of a particular \code{type} if any) within the
#' current .Rdata file.  If the call to this function (or one like it) is
#' called from within a function that runs a series of analyses then the
#' collection is limited to the function frame (i.e., only models defined
#' within the function). Thus, it is wise to either use a different .Rdata file
#' for each data set (e.g., one for dipper, another for edwards.eberhardt, etc)
#' or to run everything within functions as illustrated by \code{\link{dipper}}
#' or \code{\link{edwards.eberhardt}}. Using a separate .Rdata file is
#' equivalent to having separate .DBF/.FPT files with MARK. It is important to
#' note that functions such as \code{\link{adjust.chat}} will adjust the value
#' of chat across analyses unless specifically given a list of models to
#' adjust.
#' @export
#' @author Jeff Laake
#' @seealso \code{\link{merge.mark}},\code{\link{remove.mark}},
#' \code{\link{collect.model.names}},\code{\link{run.models}},\code{\link{model.table}},\code{\link{dipper}}
#' @keywords utility
#' @examples
#' 
#' # see examples in dipper, edwards.eberhardt and example.data
#' 
collect.models <-
function(lx=NULL,type=NULL,table=TRUE,adjust=TRUE,external=FALSE)
{
#
# collect.models - collects models in list lx of specified type (if any)
#                  and returns models in a list with a table of model
#                  results if table=TRUE
#
# Arguments:
#
#  lx       - if NULL, constructs list from frame of calling function
#             otherwise it uses list of names
#
#  type     - a character string specifying type of mark models (e.g., "CJS");
#             if NULL, all types are used
# 
#  table    - if TRUE, a table of model results is also included in list
#
#  adjust   - if TRUE, adjusts number of parameters to full rank (# of columns of
#              design matrix.
#  external - if TRUE the mark objects are saved externally rather than in the 
#              resulting marklist; the filename for each object is kept in its place
#
#  Value: 
#
#  model.list - a list of mark model objects and a table of model results
#               if table=TRUE
#
# Functions used: collect.model.names, model.table
#
#
# If lx=NULL, collect names of objects in parent frame
#
if(is.null(lx))lx=ls(envir=parent.frame())
#
# Collect names of mark models of specified type (if any) from lx
# 9 Jan 06; pulled stop from collect.model.names and put it here
#
x=collect.model.names(lx,type)
if(is.null(x))stop()
#
# Create list of mark model objects
#
z=eval(parse(text=paste("list(",paste(paste(x,"=",x,sep=""),collapse=","),")")),envir=parent.frame())
#
# Handle external objects
#
if(external)
  for (i in 1:length(z))
    if(is.list(z[[i]]))
    {
       model=z[[i]]
       marksave=paste(z[[i]]$output,".rda",sep="")
       save(model,file=marksave)                                  
       class(marksave)=class(model)
       z[[i]]=marksave
    }
#
# Return list with model table if requested
# 9 Jan 06; reordered model.table in return; also gave class of marklist to
# enable cleanup function
#
if(table)
    z$model.table=model.table(x,type,pf=2,adjust=adjust)
class(z)="marklist"
return(z)
}

