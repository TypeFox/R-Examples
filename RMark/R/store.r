#' Store models externally or restore to workspace from external storage
#' 
#' Stores/restores all mark model objects in a marklist either to or from
#' external storage.
#' 
#' For \code{store}, each mark model is stored externally and the object in the
#' list is replaced with the filename of the object. \code{restore} does the
#' opposite of storing the saved external object into the marklist and then
#' deleting the saved file.
#' 
#' @aliases store restore
#' @param x marklist of models
#' @return A modified marklist to replace the previous marklist specified as
#' the argument.
#' @author Jeff Laake
#' @export store restore
#' @keywords utility
store=function(x)
{
  if(class(x)!="marklist") stop("\nThis function only works on a marklist\n")
  for (i in 1:(length(x)-1))
  {
    model=x[[i]]
    save.mark=paste(model$output,".rda",sep="")
    save(model,file=save.mark)
    saveclass=class(x[[i]])
    x[[i]]=save.mark
    class(x[[i]])=saveclass
  }
return(x)
}
restore=function(x)
{
  model=NULL
  if(class(x)!="marklist") stop("\nThis function only works on a marklist\n")
  for (i in 1:(length(x)-1))
  {
    load(file=x[[i]])
    unlink(x[[i]])
    x[[i]]=model
  }
return(x)
}
