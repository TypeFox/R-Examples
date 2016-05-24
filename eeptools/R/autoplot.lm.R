##' A function to replicate the basic plot function for linear models in ggplot2
##' @description This uses ggplot2 to replicate the plot functionality for lm 
##' in ggplot2 and allow themes.
##' @param object a linear model object from \code{\link{lm}}
##' @param which which of the tests do we want to display output from
##' @param mfrow Describes the layout of the resulting function in the plot frames
##' @param ... additional parameters to pass through
##' @return A ggplot2 object that mimics the functionality of a plot of linear model.
##' @references Modified from: http://librestats.com/2012/06/11/autoplot-graphical-methods-with-ggplot2/
##' @seealso \code{\link{plot.lm}} which this function mimics
##' @export
##' @import ggplot2
##' @examples
##' # Univariate
##' a <- runif(1000)
##' b <- 7 * a + rnorm(1)
##' mymod <- lm(b~a)
##' autoplot(mymod)
##' # Multivariate
##' data(mpg)
##' mymod <- lm(cty~displ + cyl + drv, data=mpg)
##' autoplot(mymod)
##' 
autoplot.lm <- function(object, which=c(1:6), mfrow=c(3,2), ...){
  df <- ggplot2::fortify(object)
  df <- cbind(df, rows=1:nrow(df))
  # residuals vs fitted
  g1 <- ggplot(df, aes(.fitted, .resid)) +
    geom_point()  +
    geom_smooth(se=FALSE, method = "loess") +
    geom_hline(yintercept = 0, linetype=2, size=.2) +
    scale_x_continuous("Fitted Values") +
    scale_y_continuous("Residual") +
    labs(title="Residuals vs Fitted")+
    theme_dpi()
  # normal qq
  a <- quantile(df$.stdresid, c(0.25, 0.75))
  b <- qnorm(c(0.25, 0.75))
  slope <- diff(a)/diff(b)
  int <- a[1] - slope * b[1]
  g2 <- ggplot(df, aes(sample=.stdresid)) +
    stat_qq() +
    geom_abline(slope=slope, intercept=int) +
    scale_x_continuous("Theoretical Quantiles") +
    scale_y_continuous("Standardized Residuals") +
    labs(title="Normal Q-Q")+theme_dpi()
  # scale-location
  g3 <- ggplot(df, aes(.fitted, sqrt(abs(.stdresid)))) +
    geom_point() +
    geom_smooth(se=FALSE, method = "loess") +
    scale_x_continuous("Fitted Values") +
    scale_y_continuous("Root of Standardized Residuals") +
    labs(title="Scale-Location")+theme_dpi()
  # cook's distance
  g4 <-  ggplot(df, aes(rows, .cooksd, ymin=0, ymax=.cooksd)) +
    geom_point() + geom_linerange() +
    scale_x_continuous("Observation Number") +
    scale_y_continuous("Cook's distance") +
    labs(title="Cook's Distance")+theme_dpi()
  # residuals vs leverage
  g5 <- ggplot(df, aes(.hat, .stdresid)) +
    geom_point() +
    geom_smooth(se=FALSE, method = "loess") +
    geom_hline(yintercept = 0, linetype=2, size=.2) +
    scale_x_continuous("Leverage") +
    scale_y_continuous("Standardized Residuals") +
    labs(title="Residuals vs Leverage")+theme_dpi()
  # cooksd vs leverage
  g6 <- ggplot(df, aes(.hat, .cooksd)) +
    geom_point() +
    geom_smooth(se=FALSE, method = "loess") +
    scale_x_continuous("Leverage") +
    scale_y_continuous("Cook's distance") +
    labs(title="Cook's dist vs Leverage")+theme_dpi()
  plots <- list(g1, g2, g3, g4, g5, g6)
  # making the plots
  grid::grid.newpage()
  if (prod(mfrow)>1) {
    mypos <- expand.grid(1:mfrow[1], 1:mfrow[2])
    mypos <- mypos[with(mypos, order(Var1)), ]
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(mfrow[1], mfrow[2])))
    formatter <- function(.){}
  } else {
    mypos <- data.frame(matrix(1, length(which), 2))
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 1)))
    formatter <- function(.) {
      .dontcare <- readline("Hit <Return> to see next plot: ")
      grid::grid.newpage()
    }
  }
  j <- 1
  for (i in which){
    formatter()
    print(plots[[i]], vp=grid::viewport(layout.pos.row=mypos[j,][1], 
                                        layout.pos.col=mypos[j,][2]))
    j <- j+1
  }
}
utils::globalVariables(c(".fitted", ".resid",".stdresid",".cooksd","rows",".hat"))