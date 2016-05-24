#' Plot Construction
#' 
#' \code{"+.gg"} is a local copy of the ggplot2 add function, no change other than exporting from the ggtern namespace
#' @param e1 first object
#' @param e2 second object
#' @rdname plotconstruction
#' @export
"+.gg" <- function(e1, e2){  
  # Get the name of what was passed in as e2, and pass along so that it
  # can be displayed in error messages
  e2name <- deparse(substitute(e2))
  if      (is.theme(e1))  add_theme( e1, e2, e2name)
  else if (is.ggplot(e1)) ggint$add_ggplot(e1, e2, e2name)
}
