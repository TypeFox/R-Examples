#' Create an empty plot with compact axis notation
#' 
#' The \code{eplot()} function draws an empty plot to which the user can
#' add points, lines, text, etc. The axis notation is more compact than the
#' defaults for the \code{plot()} function. Also, axis and label are
#' appropriately suppressed when the plot occurs as part of a matrix. The
#' \code{aplot()} function simply calls \code{eplot()} again, using the
#' same arguments (with the exception of \code{main}).
#' 
#' This function simply draws an empty plot with compact axis notation, to
#' which the user can add points, lines, text, and so on. Also, if the plot
#' appears as part of a matrix, the x-axis notation is suppressed unless the
#' plot appears along the bottom row and the y-axis notation is suppress unless
#' the plot appears along the left column.
#' 
#' @aliases eplot
#' @param xlim the x limits (x1, x2) of the plot.
#' @param ylim the y limits of the plot.
#' @param xlab a label for the x axis, defaults to empty space.
#' @param ylab a label for the y axis, defaults to empty space.
#' @param main a label for the subplot. Intended for labeling a each plot in a
#' matrix. If you need a title for the entire matrix of plots, or a single
#' plot, I recommend using a call to the \code{mtext()} function.
#' @param text.size a numerical value giving the amount by which axis notation
#' should be magnified.  Reasonable values range from about 0.5 to 2.
#' @param tick.length the length of tick marks as a fraction of the smaller of
#' the width or height of the plotting region. Reasonable values range from
#' about 0.01 to 0.1.
#' @param xpos,ypos controls the distance from the tick labels to the axis.
#' Reasonable values range from about -1 to 1.
#' @param xat,yat the location of the tick marks along the axes. If "none," then
#' the axis will not be annotated.
#' @param xticklab,yticklab the labels for the tick marks. A character vector
#' the length of \code{xat} and \code{yat}.
#' @param xlabpos,ylabpos controls the distance from the axis labels to the
#' axes. Reasonable values range from about 1 to 3.
#' @param annx,anny include annotations for x and y axes?
#' @param box should a box be plotted?
#@note %% ~~further notes~~
#' @author Carlisle Rainey (\href{mailto:carlislerainey@@gmail.com}{e-mail},
#' \href{http://www.carlislerainey.com}{website})
#@seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#@references %% ~put references to the literature/web site here ~
#@keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ### Plot 0: illustrating the purpose
#' 
#' # run these lines one at a time to see what happens
#' par(mfrow = c(2,2))
#' eplot(xlim = c(-1, 1), ylim = c(0, 10))
#' aplot(main = "Hey Look! No axis labels.")
#' aplot(main = "But this one has them?!")
#' aplot(main = "And this one does just what you'd expect!")
#' # after a call to eplot() or aplot(), I just add
#' # whatever I want to the plot.
#' 
#' 
#' ### Plot 1: a simple scatter plot
#' 
#' set.seed(1234)
#' x <- rnorm(100)
#' y <- x + rnorm(100)
#' 
#' par(mfrow = c(1,1), mar = c(3,3,1,1), oma = c(0,0,2,0))
#' eplot(xlim = c(min(x), max(x)), ylim = c(min(y), max(y)),
#'            xlab = "Explanatory Variable", ylab = "Outcome Variable")
#' points(x, y)
#' abline(lm(y ~ x), lwd = 3, col = "red")
#' mtext("A Clever Title", outer = TRUE)
#' 
#' 
#' ### Plot 2: a matrix of scatter plots
#' 
#' # simulation multilevel data
#' set.seed(1234)
#' group <- rep(1:11, each = 15)
#' a <- rnorm(length(unique(group)), sd = 1)
#' b <- rnorm(length(unique(group)), mean = 1, sd = .3)
#' x <- rnorm(length(group))
#' y <- a[group] + x*b[group] + rnorm(length(group))
#' 
#' ## estimate random effects models and pull out the estimates
#' #library(lme4)
#' #hier <- lmer(y ~ x + (1 + x | group))
#' #a.hat <- fixef(hier)[1] + ranef(hier)$group[, 1]
#' #b.hat <- fixef(hier)[2] + ranef(hier)$group[, 2]
#' 
#' # draw plot
#' par(mfrow = c(3,4), mar = c(.75,.75,.75,.75), oma = c(4,4,4,1))
#' for (i in 1:11) {
#'   eplot(xlim = c(min(x), max(x)), ylim = c(min(y), max(y)),
#'              xlab = "Explanatory Variable", ylab = "Outcome Variable",
#'              main = paste("Group", i))
#'   points(x[group == i], y[group == i])
#'   #abline(a = a.hat[i], b = b.hat[i])
#'   abline(lm(y[group == i] ~ x[group == i]), lty = 3)
#' }
#' 
#' # add an overall title
#' mtext("Comparing Partial Pooling and No Pooling", outer = TRUE, line = 2)
#' 
#' 
#' ### Plot 3: a matrix of scatter plots using aplot() and addxaxis()
#' 
#' # use the same estimates as before
#' 
#' # draw the first plot with eplot()
#' par(mfrow = c(3,4), mar = c(.75,.75,.75,.75), oma = c(4,4,4,1))
#' eplot(xlim = c(min(x), max(x)), ylim = c(min(y), max(y)),
#'            xlab = "Explanatory Variable", ylab = "Outcome Variable",
#'            main = "Group 1")
#' 
#' # then add stuff
#' points(x[group == 1], y[group == 1])
#' #abline(a = a.hat[1], b = b.hat[1])
#' abline(lm(y[group == 1] ~ x[group == 1]), lty = 3)
#' legend(par("usr")[1], par("usr")[4],
#'        legend = c("partial pooling", "no pooling"), lty = c(1, 3),
#'        bty = "n", bg = NA, cex = .8)
#' 
#' # draw the rest with aplot()
#' for (i in 2:11) {
#'   aplot(main = paste("Group", i))
#'   # since we don't plan to have bottom right plot,
#'   # let's add an axis to the one above
#'   if (i == 9) { addxaxis() }
#'   points(x[group == i], y[group == i])
#'   #abline(a = a.hat[i], b = b.hat[i])
#'   abline(lm(y[group == i] ~ x[group == i]), lty = 3)
#' }
#' mtext("Comparing Partial Pooling and No Pooling", outer = TRUE, line = 2)
#' 
#' @export eplot
#' 


eplot <-
  function(xlim, ylim, xlab = NULL, ylab = NULL, 
           main = NULL, text.size = 1, tick.length = 0.02,
           xpos = -.7, ypos = -.5, xat = NULL, yat = NULL,
           xticklab = NULL, yticklab = NULL,
           xlabpos = 1.5, ylabpos = 1.5, 
           annx = TRUE, anny = TRUE, 
           box = TRUE) {
    # create an empty plot
    plot(NULL, xlim = xlim, ylim = ylim, axes = F, xlab = NA, ylab = NA)
    
    # add a box
    box()
    
    # calculate adjustment factor for axis labels if the plot is a matrix
    deflate <- 1
    if (par("mfg")[3] == 2 & 
          par("mfg")[4] == 2) {
      deflate <- 0.83
    }
    if (par("mfg")[3] > 2 | 
          par("mfg")[4] > 2) {
      deflate <- 0.66
    }
    
    # Calculate the position of axis labels.
    if (length(xat) == 0) {
      xat <- pretty(xlim, n = 5)
    }
    if (length(yat) == 0) {
      yat <- pretty(ylim, n = 4)
    }
    
    # add the x axis
    if (par("mfg")[1] == par("mfg")[3] & annx == TRUE) {
      axis(side = 1, at = xat, labels = NA, tck = -tick.length, lwd = 0, lwd.ticks = 1)
      axis(side = 1, at = xat, tick = FALSE, line = xpos, cex.axis =  .9*text.size,
           labels = xticklab)
      mtext(side = 1, xlab, line = xlabpos, cex = 1*text.size*deflate)
    }
    
    # add the y axis
    if (par("mfg")[2] == 1 & anny == TRUE) {
      axis(side = 2, at = yat, las = 1, labels = NA, tck = -tick.length, lwd = 0, lwd.ticks = 1)
      axis(side = 2, at = yat, las = 1, tick = FALSE, line = ypos, cex.axis =  .9*text.size,
           labels = yticklab)
      mtext(side = 2, ylab, line = ylabpos, cex = 1*text.size*deflate)
    }
    
    # add the plot label
    mtext(side = 3, main, line = .1, cex = 1*text.size*deflate)
    
    #   plotPar <<- list(xlim = xlim, ylim = ylim, 
    #                   xlab = xlab, ylab = ylab, 
    #                   main = main, text.size = text.size,
    #                   tick.length = tick.length, 
    #                   xpos = xpos, ypos = ypos,
    #                   xat = xat, yat = yat,
    #                   xlabpos = xlabpos, ylabpos = ylabpos,
    #                   annx = annx, anny = anny,
    #                   box = box)
    .compactrEnv$plotPar <- list(xlim = xlim, ylim = ylim, 
                                 xlab = xlab, ylab = ylab, 
                                 main = main, text.size = text.size,
                                 tick.length = tick.length, 
                                 xpos = xpos, ypos = ypos,
                                 xat = xat, yat = yat,
                                 xticklab = xticklab, yticklab = yticklab,
                                 xlabpos = xlabpos, ylabpos = ylabpos,
                                 annx = annx, anny = anny,
                                 box = box)
  }


.compactrEnv <- new.env()
