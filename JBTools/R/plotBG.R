plotBG <- function(
  ##title<< Plot a colored plot background
  color         ##<< color: which color to use
  ,exp.factor=1 ##<< integer:possible expansion factor to use if very high values are plotted
  , xlim = c()
  , ylim = c()
  ,...          ##<< further arguments passed to plot()
  )
  ##description<< plotBG colors the plot background of the plot region (not the device region!).
  ##details<<
  ## The function opens a plot with no content and plots a very large polygon of the given
  ## color. Afterwards a new high level plot can be plotted over this background.
{
  if (length(xlim) == 0)
    xlim = exp.factor*c(-1e20, 1e20)
  if (length(ylim) == 0)
    ylim = exp.factor*c(-1e20, 1e20)
  plot(1, 1, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', type = 'n',
       xlim = xlim, ylim = ylim, ...)
  polygon(c(par()$usr[1:2], rev(par()$usr[1:2])), rep(par()$usr[3:4], each = 2), col = color)
  par(new=TRUE)
   ##value<< Nothing is returned.
}
