#' @title Exploratory Graphs for Two Factor Designs
#' 
#' @description Function creates side-by-side boxplots for each factor, a design plot (means), and an interaction plot.
#' 
#' @param Y response variable 
#' @param fac1 factor one
#' @param fac2 factor two
#' @param COL a vector with two colors
#' 
#' @author Alan T. Arnholt <arnholtat@@appstate.edu>
#' 
#' @seealso \code{\link{oneway.plots}}, \code{\link{checking.plots}}
#' 
#' @export
#' 
#' @examples
#' with(data = TIREWEAR, twoway.plots(wear, treat, block))
#' #################################
#' ## Similar graphs with ggplot2 ##
#' #################################
#' p1 <- ggplot(data = TIREWEAR, aes(x = treat, y = wear, fill = treat)) +
#' geom_boxplot() + guides(fill = FALSE) + theme_bw()
#' p2 <- ggplot(data = TIREWEAR, aes(x = block, y = wear, fill = block)) +
#' geom_boxplot() + guides(fill = FALSE) + theme_bw() 
#' p3 <- ggplot(data = TIREWEAR, aes(x = treat, y = wear, color = block,
#' group = block)) + stat_summary(fun.y = mean, geom = "point", size = 4) + 
#' stat_summary(fun.y = mean, geom = "line") + theme_bw()
#' p4 <- ggplot(data = TIREWEAR, aes(x = treat, y = wear, color = treat)) + 
#' geom_boxplot() + facet_grid(. ~ block) +theme_bw()
#' p1
#' p2
#' p3
#' p4
#' ## To get all plots on the same device use gridExtra (not run)
#' ## library(gridExtra)
#' ## grid.arrange(p1, p2, p3, p4, nrow=2)
#' @keywords hplot
####################################################################
twoway.plots<-function(Y, fac1, fac2, COL=c("#A9E2FF", "#0080FF")){
  opar <- par(no.readonly = TRUE)
  par(mfrow=c(2, 2), mar = c(5.1, 4.1, 1.1, 1.1))
  YL <- range(Y)
  plot(Y ~ fac1, col = COL[1], xlab = deparse(substitute(fac1)),
       ylab = deparse(substitute(Y)), ylim = YL)
  plot(Y ~ fac2, col = COL[2], xlab = deparse(substitute(fac2)),
       ylab = deparse(substitute(Y)), ylim = YL)
  plot.design(Y ~ fac1 + fac2, fun = "mean", 
              ylab = deparse(substitute(Y)), ylim = YL)
  interaction.plot(fac1, fac2, Y, xlab = deparse(substitute(fac1)), 
                   trace.label = deparse(substitute(fac2)),
                   type = "b", legend = FALSE, 
                   ylab = deparse(substitute(Y)), ylim = YL)
  on.exit(par(opar))
}
