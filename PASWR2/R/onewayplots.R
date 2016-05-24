#' @title Exploratory Graphs for Single Factor Designs
#' 
#' @description Function to create dotplots, boxplots, and design plot (means) for single factor designs
#' 
#' @param Y response variable for a single factor design
#' @param fac1 predictor variable (factor)
#' @param COL a vector with two colors
#' 
#' @author Alan T. Arnholt <arnholtat@@appstate.edu>
#' 
#' @seealso \code{\link{twoway.plots}}, \code{\link{checking.plots}}
#' 
#' @export
#' 
#' @examples
#' with(data = TIRE, oneway.plots(stopdist, tire))
#' ## Similar graphs with ggplot2 
#' ggplot(data = TIRE, aes(tire, stopdist, fill = tire)) + 
#' geom_dotplot(binaxis = "y", stackdir = "center") + coord_flip() + theme_bw()
#' ggplot(data = TIRE, aes(tire, stopdist, fill = tire)) + geom_boxplot() +
#' guides(fill = FALSE) + theme_bw()
#' 
#' @keywords hplot
#############################################################
oneway.plots <- function(Y, fac1, COL = c("#A9E2FF","#0080FF")){
  opar <- par(no.readonly = TRUE)
  par(mfrow = c(2, 2), pty = "m",mar = c(4.1, 2.1, 1.1, 1.1))
  YL <- range(Y)
  dotchart(x = Y, groups = fac1, color = COL[2], pch = 1,
           xlab = deparse(substitute(Y)), xlim = YL,
           gdata = tapply(Y, fac1, mean), gpch = 17)
  boxplot(Y, col = COL[1], xlab = deparse(substitute(Y)), horizontal = TRUE)
  boxplot(Y ~ fac1, col = COL[1], horizontal = FALSE, xlab = deparse(substitute(fac1)))
  plot.design(Y ~ fac1,  ylim = YL, xlab = "Main Factor")
  on.exit(par(opar))
}
