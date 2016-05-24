
#' @import ggplot2
NULL

#' Example grid
#'
#' An example grid object made of modified Iris-data with package preprocomb.
#'
#' modifiediris <- droplevels(iris[-c(1:60),])
#'
#' examplegrid <- preprocomb::setgrid(phases=c("scaling", "smoothing", "outliers", "selection", "sampling"), data=modifiediris)
#' @format A GridClass object
"examplegrid"

getuniquepreprocessors <- function(x) {factor(unlist(unique(x)))}

globalVariables(c("value","Var1", "examplegrid"))


