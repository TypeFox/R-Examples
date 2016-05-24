# Need to add labels as option


#' Display PIM for a parameter
#' 
#' Extract PIMS for a particular parameter and display either the full PIM
#' structure or the simplified PIM structure.
#' 
#' 
#' @param model mark model object
#' @param parameter character string of a particular type of parameter in the
#' model (eg "p","Phi","pent","S")
#' @param simplified if TRUE show simplified PIM structure; otherwise show full
#' structure
#' @param use.labels if TRUE, uses time and cohort labels for columns and rows
#' respectively
#' @return None
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{make.design.data}}
#' @keywords utility
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(dipper)
#' results=mark(dipper)
#' PIMS(results,"Phi")
#' PIMS(results,"Phi",simplified=FALSE)
#' } 
PIMS=function(model,parameter,simplified=TRUE,use.labels=TRUE)
{
#
# Extract PIMS for a particular parameter and display either the full PIM structure
# or the simplified PIM structure
#
# Arguments:
#
#  model     - mark model object
#  parameter - character string of a particular type of parameter in the model (eg "p","Phi","pent","S")
#  simplified- if TRUE show simplified PIM structure; otherwise show full structure
#
  model=load.model(model)
  if(!valid.parameters(model$model,parameter)) stop()
   number.of.pims = length(model$pims[[parameter]])
   nn = names(model$pims[[parameter]][[1]])
   nn = nn[nn != "pim"]
   for (i in 1:number.of.pims) {
       struct = model$pims[[parameter]][[i]][2:(length(nn) +
           1)]
       struct[nn == "group"] = model$group.labels[unlist(struct[nn ==
           "group"])]
       cat(paste(nn, "=", struct, collapse = ";"), "\n")
       pim = model$pims[[parameter]][[i]]$pim
       fullpim = pim
       if (simplified)
           pim[pim > 0] = model$simplify$pim.translation[pim[pim >
               0]]
       pim[pim == 0] = NA
       if (!use.labels) {
           row.names(pim) = 1:dim(pim)[1]
           colnames(pim) = 1:dim(pim)[2]
       }
       else {
           colnames(pim) = model$design.data[[parameter]]$time[1:dim(pim)[2]]
           if (!is.null(model$design.data[[parameter]]$cohort))
               row.names(pim) = model$design.data[[parameter]]$cohort[(diag(fullpim) -
                 fullpim[1, 1] + 1)]
           else row.names(pim) = 1:dim(pim)[1]
       }
       print.table(pim, na.print = "")
   }
}
