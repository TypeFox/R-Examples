#' Print MARK objects
#' 
#' Displays MARK output file or input file
#' with MarkViewer (notepad.exe by default) so it can be viewed. 
#' 
#' If the model has been run (\code{model$output} exists) the output file
#' stored in the directory as identified by the basefile name
#' (\code{model$output}) and the suffix ".out" is displayed with a call to
#' \code{MarkViewer}. If \code{input} is set to \code{TRUE} then the MARK input
#' file is displayed instead. By default the \code{MarkViewer} is notepad but
#' any program can be used in its place that accepts the filename as the first
#' argument. For example setting \code{MarkViewer="wp"} will use wordperfect
#' (wp.exe) as long as wp.exe is in the search path.  \code{MarkViewer} must be
#' set during each R session, so it is best to include it in your \code{.First}
#' function to change it permanently.  Since \code{print.mark} is the generic
#' function to print \code{mark} objects you can use it by just typing the name
#' of a \code{mark} object at the R prompt and it will call \code{print.mark}.
#' For example, if \code{mod} is a \code{mark} object then typing \code{mod} is
#' the same as \code{print.mark(mod)}
#' 
#' @usage \method{print}{mark}(x,...,input=FALSE)
#' @param x mark model object; or list of mark model objects created with
#' \code{\link{collect.models}}
#' @param ... additional non-specified argument for S3 generic function
#' @param input if TRUE, prints mark input file; otherwise the output file
#' @return None
#' @export 
#' @author Jeff Laake
#' @seealso \code{\link{summary.mark}}
#' @keywords utility
print.mark <- function(x,...,input=FALSE)
{
# -------------------------------------------------------------------------------------------------------------
#
# print.mark  - extracts model output file from MARK and sends to notepad so it can be viewed
#
# Arguments:
#
#  x - a mark model object
#  input  - if TRUE, prints the mark input file rather than the output file
#
# Value:
#  None
#
# -------------------------------------------------------------------------------------------------------------
#
#
#
   os=R.Version()$os
   if(!exists("MarkViewer"))
     if(os=="mingw32")
        MarkViewer="notepad"
     else
        MarkViewer="pico"
#
#  If the model has been run (model$output exists) extract to temp file dummy.xxx and 
#  call notepad to view it; any program could be used in its place. The temp file is immediately deleted
#  9 Jan 06; changed such that model$output is just the baseline value of the filename
#
#  def.options=options()
#  options(useFancyQuotes=FALSE)
  model=load.model(x)
  if(!input)
  {
      if(!is.null(model$output))
      {
         if(file.exists(paste(model$output,".out",sep="")))
         {
            if(os=="mingw32")
               system(paste(shQuote(MarkViewer),paste(model$output,".out",sep="")),invisible=FALSE,wait=FALSE)
            else
               system(paste(MarkViewer,paste(model$output,".out",sep="")),wait=FALSE)
         }
         else
            cat(paste("Cannot locate file ",model$output,".out\n",sep=""))
      }else
        print.default(model)
  }
  else
  {
      if(!is.null(model$output))
      {
         if(file.exists(paste(model$output,".inp",sep="")))
         {
            if(os=="mingw32")
               system(paste(shQuote(MarkViewer),paste(model$output,".inp",sep="")),invisible=FALSE,wait=FALSE)
            else
               system(paste(MarkViewer,paste(model$output,".inp",sep="")),wait=FALSE)
         }
         else
            cat(paste("Cannot locate file ",model$output,".inp\n",sep=""))
      }else
        print.default(model)
  }

#   options(def.options)
   invisible()
}
#' Print MARK list objects
#' 
#' Displays the model.table if it exists. To display the
#' output for a mark model contained in a list, simply type the list value
#' (e.g., typing mymarklist[[2]] will display output for the second model).
#' The function \code{print.marklist} was created to avoid accidental typing of
#' the model list which would call \code{print.mark} for each of the models.
#' 
#' If the model has been run (\code{model$output} exists) the output file
#' stored in the directory as identified by the basefile name
#' (\code{model$output}) and the suffix ".out" is displayed with a call to
#' \code{MarkViewer}. If \code{input} is set to \code{TRUE} then the MARK input
#' file is displayed instead. By default the \code{MarkViewer} is notepad but
#' any program can be used in its place that accepts the filename as the first
#' argument. For example setting \code{MarkViewer="wp"} will use wordperfect
#' (wp.exe) as long as wp.exe is in the search path.  \code{MarkViewer} must be
#' set during each R session, so it is best to include it in your \code{.First}
#' function to change it permanently.  Since \code{print.mark} is the generic
#' function to print \code{mark} objects you can use it by just typing the name
#' of a \code{mark} object at the R prompt and it will call \code{print.mark}.
#' For example, if \code{mod} is a \code{mark} object then typing \code{mod} is
#' the same as \code{print.mark(mod)}
#' 
#' @aliases print.marklist
#' @usage  \method{print}{marklist}(x,...)
#' @param x mark model object; or list of mark model objects created with
#' \code{\link{collect.models}}
#' @param ... additional non-specified argument for S3 generic function
#' @return None
#' @export 
#' @export print.marklist
#' @author Jeff Laake
#' @seealso \code{\link{summary.mark}}
#' @keywords utility

print.marklist<-function(x,...)
{
   ncol=dim(x$model.table)[2]
   if(!is.null(x$model.table))
   {
     if(is.null(x$model.table$chat))
        print(x$model.table[,(ncol-5):ncol])
     else
        print(x$model.table[,(ncol-6):ncol])
   } else cat("No model.table is available")
}

