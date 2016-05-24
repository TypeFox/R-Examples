#' bartab
#'
#' Plot the overlap of three groups with a barplot
#'
#' @param x logical
#' @param y logical
#' @param z logical
#' @param names a character vector of length 3
#' @param skipNone remove the "none" group
#' @param ... further arguments passed on to \code{\link{barplot}}
#' @author Michael I. Love
#' 
#' @examples
#' 
#' set.seed(1)
#' x <- sample(c(FALSE,TRUE), 10, replace=TRUE)
#' y <- sample(c(FALSE,TRUE), 10, replace=TRUE)
#' z <- sample(c(FALSE,TRUE), 10, replace=TRUE)
#' bartab(x,y,z,c("X","Y","Z"))
#' 
#' 
bartab <- function(x,y,z,names,skipNone=FALSE,...) {
  x <- factor(x,c("FALSE","TRUE"))
  y <- factor(y,c("FALSE","TRUE"))
  z <- factor(z,c("FALSE","TRUE"))
  tabs <- as.vector(table(x, y, z))
  names(tabs) <- c("none",names[1],names[2],paste(names[1:2],collapse="+"),
                   names[3],paste(names[c(1,3)],collapse="+"),
                   paste(names[2:3],collapse="+"),"all")
  tabs <- tabs[c(1,2,3,5,4,6,7,8)]
  if (skipNone) {
    tabs <- tabs[-1]
  }
  barplot(tabs, ...)
}
