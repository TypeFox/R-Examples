#' A palette of 17 diverging colors
#' 
#' 17 diverging colors created by combining the Set1 and Dark2 palettes from 
#' RColorBrewer.
#' 
#' @source R Color brewer: \code{c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"))}.
#' @format Vector of 17 diverging colors.
#' @examples
#' \dontrun{
#'  barplot(rep(1, 17), col=Tcols, axes=F, names=c(rep("Set1", 9), rep("Dark2", 8)), horiz=TRUE, las=2)
#' }

"Tcols"