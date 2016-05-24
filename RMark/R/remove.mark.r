#' Remove mark models from list
#' 
#' Remove one or more mark models from a marklist
#' 
#' 
#' @param marklist an object of class "marklist" created by
#' \code{\link{collect.models}} or \code{\link{merge.mark}}
#' @param model.numbers vector of one more model numbers to remove from the
#' marklist
#' @return model.list: a list of \code{mark} models and a table of model
#' results.
#' @author Jeff Laake
#' @export
#' @seealso
#' \code{\link{collect.models}},\code{\link{merge.mark}},\code{\link{run.models}},\code{\link{model.table}},\code{\link{dipper}}
#' @keywords utility
#' @examples
#' 
#' # see example in dipper
#' 
remove.mark=function(marklist,model.numbers)
{
#
# Removes a set of models from a marklist and returns the modified
# marklist
#
# Check validity of arguments
#
  if(class(marklist)[1]!="marklist")
     stop(paste("\n",substitute(marklist), "is not a marklist object\n"))
  if(!is.null(marklist$model.table))
     marklist$model.table=NULL
  if(max(model.numbers)>length(marklist) | min(model.numbers)<1)
     stop("\nInvalid set of model numbers for removal\n")
#
# Remove specified models
#
  marklist=marklist[-model.numbers]
#
# If any models are left, create a model table otherwise return NULL
#
  if(length(marklist)>0)
  {
     marklist$model.table=model.table(marklist)
     class(marklist)="marklist"
     return(marklist)
  }
  else
     return(NULL)
}
