
#' Plot the history of the scores of each team over time.
#' 
#' The best score of each team has a bold symbol.
#' 
#' @param history   list of the submissions history per team as returned by \code{\link{compute_metrics}}
#' @param metric    string. name of the metric considered
#' @param test_name string. name of the test set used: \code{"quiz"} or \code{"test"}
#' @param baseline  string. name of the team considered as the baseline. Its best score
#'   will be plotted as a constant and will not appear in the legend.
#' @param col       colors of the teams
#' @param pch       symbols of the teams
#' @param by        real. interval width of grid lines
#' @param xlab,ylab axis labels. see \code{\link[graphics]{title}}.
#' @param ...       further parameters passed to \code{\link[graphics]{plot}} function.
#' @param bty,fg,col.axis,col.lab graphical parameters. see \code{\link[graphics]{par}}.
#' @param text.col the color used for the legend text. see \code{\link[graphics]{legend}}.
#' 
#' @importFrom graphics plot axis.POSIXct axTicks abline lines points legend
#' @importFrom grDevices palette adjustcolor
#' @importFrom stats dnorm density
#' @export
#' @return \code{NULL}
plot_history <- function(history, metric, test_name="quiz", baseline="baseline", 
                         col=1:length(history), pch=rep(21:25, 100), by = .05,
                         xlab="Date", ylab = "Score", bty='l',
                         fg="darkslategray", col.axis=fg, col.lab=fg, 
                         text.col = fg, ...) {
  metric_column = paste(metric, test_name, sep=".")
  # get baseline score
  ind_base = which(names(history)==baseline)
  if (length(ind_base)>0) {
    base_score = min(history[[ind_base]][[metric_column]])
    history = history[-ind_base]
  }
  else
    base_score = NULL
  
  # compute axis limits
  if (length(history)==0)
    xlim = c(Sys.time()-1, Sys.time())
  else 
    xlim = range(history[[1]]$date)
  ylim = c(base_score-by, base_score+by)
  for (i in seq(along=history)) {
    xlim = range(xlim, history[[i]]$date)
    ylim = range(ylim, history[[i]][[metric_column]])
  }
  
  # empty figure
  plot(NA, type='n', xlab=xlab, ylab = ylab, bty=bty,
       xlim=xlim, ylim=ylim, xaxt = 'n', fg=fg, col.axis=col.axis, col.lab=col.lab, ...)
  axis.POSIXct(1, at = seq(xlim[1], xlim[2], by=3600*24*7), format="%d %b", fg=fg, col.axis=col.axis)
  
  # grid
  aty = range(axTicks(2))
  abline(h=seq(aty[1]-by, aty[2]+by, by=by), col="lightgray", lty=3)
  
  # baseline
  abline(h=base_score, col=fg, lty=2, lwd=2)
  
  # history for each team
  for (i in seq(along=history)) {
    lines(history[[i]]$date, history[[i]][[metric_column]], 
          col=col[i], lty=3)
    points(history[[i]]$date, history[[i]][[metric_column]], 
           col=col[i], pch=pch[i], lwd=2)
    # best submission in bold
    ind = which.min(history[[i]][[metric_column]])
    points(history[[i]]$date[ind], history[[i]][[metric_column]][ind], 
           col=col[i], pch=pch[i], bg=col[i], lwd=2)
  }
  
  # legend
  leg = c(baseline, names(history))
  legend('topright', legend = leg, col=c(fg,col), pch=c(NA,pch), lwd=c(2,rep(2,length(history))), 
         lty=c(2,rep(NA,length(history))), bty='n', xpd = NA, inset = c(-0.22, 0),
         text.col = text.col)
  
  invisible(NULL)
}


#' Plot the density of submissions over time.
#' 
#' @param history   list of the submissions history per team as returned by \code{\link{compute_metrics}}
#' @param baseline  string. name of the team considered as the baseline that will not be plotted.
#' @param col       colors of the teams.
#' @param alpha.f   factor modifying the opacity alpha of colors; typically in [0,1].
#' @param bw        real. the smoothing bandwidth to be used by \code{\link[stats]{density}} in seconds.
#' @param by        real. height of the interval between two teams in nb of submissions.
#' @param xlab,ylab axis labels. see \code{\link[graphics]{title}}.
#' @param ...       further parameters passed to \code{\link[graphics]{plot}} function.
#' @param bty,fg,col.axis,col.lab graphical parameters. see \code{\link[graphics]{par}}.
#' @param text.col the color used for the legend text. see \code{\link[graphics]{legend}}.
#' 
#' @return \code{NULL}
#' 
#' @importFrom graphics plot axis.POSIXct polygon legend
#' @importFrom grDevices palette adjustcolor
#' @importFrom stats dnorm density
#' @export
#' @seealso \code{\link[stats]{density}}
plot_activity <- function(history, baseline="baseline", col=1:length(history), 
                          alpha.f = .7, bw = 3600*24, by = 4,
                          xlab="Date", ylab = "Submissions density", bty='l',
                          fg="darkslategray", col.axis=fg, col.lab=fg, 
                          text.col = fg, ...) {
  # baseline index
  ind_base = which(names(history)==baseline)
  history <- history[-ind_base]
  
  # compute axis limits
  if (length(history)==0)
    xlim = c(Sys.time()-1, Sys.time())
  else 
    xlim = range(history[[1]]$date)
  
  ylim <- c(0, by*(length(history)+1))
  for (i in seq(along=history)) {
    xlim = range(xlim, history[[i]]$date)
  }
  
  # empty figure
  plot(NA, type='n', xlab=xlab, ylab = ylab, bty=bty,
       xlim=xlim, ylim=ylim, xaxt = 'n', yaxt='n', fg=fg, col.axis=col.axis, col.lab=col.lab, ...)
  axis.POSIXct(1, at = seq(xlim[1], xlim[2], by=3600*24*7), format="%d %b", fg=fg, col.axis=col.axis)
  
  # palette with transparency
  mycols <- palette()
  mycols_alpha <- adjustcolor(mycols, alpha.f)
  palette(mycols_alpha)
  
  # submissions density for each team
  ymax <- dnorm(0, sd = bw)
  for (i in seq(along=history)) {
    dens <- density(as.numeric(history[[i]]$date), bw = bw)
    dens$y <- dens$y * dens$n / ymax # rescale
    ybase <- (length(history)-i)*by
    polygon(c(dens$x, rev(dens$x)), 
            c(rep(ybase, length(dens$y)), ybase+rev(dens$y)),
            col = col[i], border=NA)
  }
  
  # legend
  if (length(history)>0)
    legend('topright', legend = names(history), col = col, pch = 19, lwd = NA, pt.cex = 1.5,
           bty='n', xpd = NA, inset = c(-0.22, 0), text.col = text.col)
  
  # restore palette
  palette(mycols)
  
  invisible(NULL)
}
