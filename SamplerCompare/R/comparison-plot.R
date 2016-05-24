# From SamplerCompare, (c) 2010 Madeleine Thompson

# comparison.plot generates a plot comparing MCMC performance as
# described in "Graphical Comparison of MCMC Performance".  See
# ?comparison.plot for usage details.

globalVariables(c('tuning', 'evals', 'act', 'act.025', 'act.975', 'lim'))

comparison.plot <- function(RS, xlab=NULL, ylab=NULL, base_size=10, ...) {
  # First, plot the results with finite ACT.

  RSfinite <- subset(RS, is.finite(RS$act))

  # Compute a reasonably-spaced set of tick marks.

  x.breaks <- log.breaks(RSfinite$tuning, 10)

  # Manually handle xlab and ylab so they can be overridden by caller.

  if (is.null(xlab))
    xlab <- 'scale tuning parameter'
  if (is.null(ylab))
    ylab <- paste('# of evals. of log density function per',
                  'uncorrelated obs. (with 95% CI)')

  # Generate grid of plots of evals*act vs. tuning parameter, with
  # a 95% confidence interval.  The theme and x scale overrides are
  # largely personal preference; if the caller wants different choices,
  # they can specify them with the + operator on the returned value.
  # "labeller" really is spelled this way in facet_grid.

  p <- ggplot2::qplot(tuning, evals*act, ymin=evals*act.025, ymax=evals*act.975,
      data=RSfinite, log='xy', geom='pointrange', xlab=xlab, ylab=ylab, ...) +
    ggplot2::facet_grid(dist.expr~sampler.expr, labeller=ggplot2::label_parsed) +
    ggplot2::scale_x_log10(breaks=x.breaks, labels=x.breaks) +
    ggplot2::theme_bw(base_size=base_size) +
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_text(angle=45, vjust=1))

  # Next, plot the results with infinite/unknown ACT as question
  # marks at lim, computed to be at the top of the plot.
  
  if (any(!is.finite(RS$act))) {
    RSinf <- subset(RS, !is.finite(RS$act))
    RSinf$lim <- max(RSfinite$evals*RSfinite$act)
    p <- p + ggplot2::geom_text(ggplot2::aes(x=tuning, y=lim, label='?'),
                                data=RSinf, size=4.0)
  }

  return(p)
}

# Return an sequence of tick marks to be used on a log axis.  Used
# to keep the automatic tick generation from spacing ticks too
# closely; this only returns ticks on integer powers of base.

log.breaks <- function(data, base) {
  breaks <- base^seq(floor(min(log(data,base=base))),
                     ceiling(max(log(data,base=base))))
}
