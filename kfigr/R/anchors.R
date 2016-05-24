#' Anchor Information
#'
#' Retrieves diagnostics such as the anchor index and history. Use for 
#' code verification and troubleshooting. Also used internally by kfigr.
#'
#' @param tag Optional specification of return type. \code{tag = "index"} 
#'   returns a dataframe listing the chunk labels, asigned types and index 
#'   numbers. \code{tag = "history"} returns a dataframe listing all reference 
#'   calls, in order.
#' @return If \code{tag = "index"}, a dataframe listing all anchored chunks. If 
#'   \code{tag = "history"}, a dataframe listing all references made, in order.
#'   If the value of \code{tag} matches a specific \code{type}, all references 
#'   of that \code{type} are provided in a list. If \code{tag} is missing, all 
#'   references of all types are provided in a nested list.
#'
#' @examples
#' figr("foo", type="figure")
#' figr("bar", type="table")
#' figr("test", type="figure")
#' anchors()
#'
#' @importFrom stats setNames
#' @export
anchors <- function(tag){
  formathist <- function(x){
    if(length(x) > 0)
      setNames(do.call(rbind.data.frame, x),
               c('label', 'type', 'number', 'referenced.by'))
    else
      x
  }
  formatindex <- function(x){
    d <- as.data.frame(cbind(names(x), x))
    rownames(d) <- NULL
    names(d) <- c("label", "type")
    d["number"] <- rep(NA, length(x))
    for(i in unique(d$type))
      d[d$type == i, "number"] = seq(1, sum(d$type == i), 1)
    d
  }
  
  types <- get('types', envir = anchorenv)
  a <- structure(setNames(lapply(types, function(x) get(x, envir=anchorenv)), types))  
  if(missing(tag))
    return(a)
  if(tag == "index")
    return(formatindex(get('index', envir=anchorenv)))
  if(tag == "history")
    return(formathist(get("history", envir=anchorenv)))
  return(a[[tag]])
}
