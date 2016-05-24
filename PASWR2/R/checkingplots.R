#' @title Checking Plots
#' 
#' @description Function that creates four graphs that can be used to help assess independence, normality, and constant variance
#' 
#' @param model an aov or lm object
#' @param n.id the number of points to identify
#' @param COL vector of two colors
#' 
#' @author Alan T. Arnholt <arnholtat@@appstate.edu>
#' 
#' @seealso \code{\link{twoway.plots}}, \code{\link{oneway.plots}}
#' 
#' 
#' @export
#' 
#' @examples
#' mod.aov <- aov(stopdist ~ tire, data = TIRE)
#' checking.plots(mod.aov)
#' rm(mod.aov)
#' 
#' # Similar graphs using ggplot2
#' #
#' mod.aov <- aov(stopdist ~ tire, data = TIRE)
#' fortify(mod.aov)
#' # library(gridExtra) used to place all graphs on the same device
#' p1 <- ggplot(data = mod.aov, aes(x = 1:dim(fortify(mod.aov))[1], y = .stdresid, 
#' color = tire)) + geom_point() + labs(y = "Standardized Residuals", 
#' x = "Ordered Residuals") + geom_hline(yintercept = c(-3,-2, 2, 3), 
#' linetype = "dashed", col = "pink") + theme_bw()
#' p2 <- ggplot(data = mod.aov, aes(sample = .stdresid, color = tire)) +
#' stat_qq()  + geom_abline(intercept = 0, slope = 1, linetype = "dashed", col = "pink") + theme_bw()
#' p3 <- ggplot(data = mod.aov, aes(x = .fitted, y = .stdresid, color = tire)) + 
#' geom_point() + geom_hline(yintercept = 0, linetype = "dashed") + 
#' labs(y = "Standardized Residuals", x = "Fitted Values") + 
#' geom_hline(yintercept = c(-3, -2, 2, 3), linetype = "dashed", color = "pink") + 
#' theme_bw()
#' p1
#' p2
#' p3
#' multiplot(p1, p2, p3, cols = 1)
#' # Or use the following (not run) to get all graphs on the same device
#' # library(gridExtra)
#' # grid.arrange(p1, p2, p3, nrow=3)
#' rm(mod.aov, p1, p2, p3)
#' 
#' @keywords hplot
#' 
############## added density plot and 2*2 design 6/27/13 ###
## Change type from "b" to "o" 11/04/14
checking.plots <- function(model, n.id = 3, COL=c("#0080FF", "#A9E2FF")){
  Identify <- function(x, y, n.id = 3, AD = 0.2){
    y <- y[!is.na(y)]
    y <- y[is.finite(y)]
    n <- length(y)
    pos <- as.numeric(names(y))
    oy <- order(abs(y))
    b <- y*AD
    W <- oy[(n - n.id + 1):n]
    text(x[W], y[W]+b[W], as.character(pos[W]))
  }
  par(mfrow=c(2, 2), pty = "m")
  varmodel <- deparse(substitute(model))
  y <- rstandard(model)
  mAx <- max(abs(y[is.finite(y)])) + .5
  n <- length(y)
  plot( (1:n), y, ylab = "standardized residuals", col = COL[1],
        xlab = "ordered values", ylim=c(-mAx, mAx), type = "o",
        main = paste("Standardized residuals versus \n ordered values for", varmodel))
  Identify((1:n), y, n.id)
  abline(h = 0, lty = 2, col = COL[2])
  abline(h = c(2, 3), col = COL, lty = 2)
  abline(h = c(-2, -3), col = COL, lty = 2)
  qqnorm(y, col = COL[1], ylab = "standardized residuals", xlim = c(-mAx, mAx),
         ylim = c(-mAx, mAx), main = paste("Normal Q-Q plot of standardized \n residuals from", varmodel))
  abline(a = 0, b = 1, col = COL[2])
  junk <- qqnorm(y, plot.it = FALSE)
  Identify(junk$x, junk$y, n.id)
  plot(fitted(model), y, col = COL[1], xlab = "fitted values", ylab = "standardized residuals",
       ylim = c(-mAx, mAx), main = paste("Standardized residuals versus \n fitted values for", varmodel))
  abline(h = 0, lty = 2, col = COL[2])
  abline(h = c(2, 3), col = COL, lty = 2)
  abline(h = c(-2, -3), col = COL, lty = 2)
  Identify(fitted(model), y, n.id)
  #
  plot(density(y), main = paste("Density plot of standardized\nresiduals for", varmodel), col = COL[1])  
  par(mfrow = c(1, 1), pty = "m")
  #  
}