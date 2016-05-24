##' Treatment Prevalence vs Duration graphic
##'
##' Provides a special graphic for treatment prevalence-vs-duration.
##' Used in place of lattice function 'panel.bwplot', this function
##' enables a special kind of graphic displaying prevalence of various
##' kinds of treatment within a patient population, together with the
##' distribution of treatment durations.
##' 
##' @param x Vector of (numeric or difftime) treatment durations
##' @param y The treatment factor
##' @param prevalence Treatment prevalence, typically an \code{xtabs}
##' object calculated by, e.g., xtabs(~trt, data)
##' @param box.width If this is specified, a more-or-less standard
##' boxplot will be produced
##' @param pch Plotting symbol used for boxplot median
##' @param pch.out Plotting symbol for outliers (default is dot)
##' @param alpha Boxplot transparency
##' @param \dots Currently used to pass \code{key.colors} parameter
##' @param stats Function to generate boxplot statistics (default is
##' \code{boxplot.stats})
##' @param coef How far boxplot whiskers should extend outside the box
##' @param do.out Logical; should outliers be plotted?
##' @author David C. Norris
##' @seealso \code{\link[lattice:panel.bwplot]{panel.bwplot}}
##' @keywords internal hplot
##' 
##' ## TODO: Document several kinds of usage, including 'passthru' to a
##' ##       standard bwplot.
##' @export panel.prevalence
panel.prevalence <-
  function (x, y, prevalence,
            box.width = NA, # May be specified in bwplot call if a standard boxplot is desired
            pch = "+",    # Plotting symbol for boxplot median
            pch.out = 20, # Plotting symbol for outliers (20 = filled dot)
            alpha = 0.5,
            ...,
            stats = boxplot.stats,
            coef = 1.5,
            do.out = TRUE) 
{
  if (is(x,"difftime"))
     x <- as.double(x, units=units(x))
  if (!is.factor(y))
    stop("Function 'panel.prevalence' requires that y be a factor.")
  ## This next requirement is a temporizing measure pending development of a sound metadata facility:
  if (is.null(list(...)$key.colors))
    stop("Function 'panel.prevalence' requires that 'key.colors' be provided as arg to bwplot.")
  if (all(is.na(x) | is.na(y))) 
    return()
  if(!is.na(box.width)){ # plot a more-or-less standard boxplot
    blist <- tapply(x, y, stats, coef = coef, do.out = do.out)
    blist <- blist[!sapply(blist, is.null)]
    blist.stats <- t(sapply(blist, "[[", "stats"))
    blist.n <- lapply(blist, "[[", "n")
    blist.out <- lapply(blist, "[[", "out")
    levels.fos <- prevalence[names(blist.n)[blist.n>0]]
    ## We must remove any degeneracy present in levels.fos,
    ## to permit their (ab)use as (real-valued!) factor levels.
    if(length(levels.fos) > length(unique(levels.fos)))
      levels.fos <- jitter(levels.fos, 0.001)
    x <- x[y %in% names(levels.fos)]
    y <- y[y %in% names(levels.fos)]
    y <- factor(as.character(y), levels=names(levels.fos))
    fill <- unlist(list(...)$key.colors)[names(levels.fos)]
    box.rectangle <- trellis.par.get("box.rectangle")
    box.umbrella <- trellis.par.get("box.umbrella")
    plot.symbol <- trellis.par.get("plot.symbol")
    trellis.par.set(list(box.rectangle = list(
                           col=fill,
                           alpha=alpha)
                         ,box.umbrella = list(
                           col=rep(fill, times=2),
                           alpha=alpha)
                         ,plot.symbol = list(
                           pch=pch.out,
                           col=rep(fill, times=sapply(blist.out[names(levels.fos)], length)),
                           alpha=alpha)
                         ))
    panel.bwplot(x, levels.fos[y], box.width = box.width,
                 pch = pch,
                 fill = fill,
                 levels.fos = as.numeric(levels.fos),
                 ...)
    trellis.par.set("box.rectangle", box.rectangle)
    trellis.par.set("box.umbrella", box.umbrella)
    trellis.par.set("plot.symbol", plot.symbol)
  } else { # simply plot all the values
    plot.symbol.restore <- trellis.par.get("plot.symbol")
    trellis.par.set(list(plot.symbol=list(pch=pch.out)))
    plot.symbol <- trellis.par.get("plot.symbol")
    panel.points(x = x,
                 y = prevalence[y],
                 pch = plot.symbol$pch,
                 col = unlist(list(...)$key.colors)[as.character(y)],
                 alpha = alpha,
                 cex = plot.symbol$cex,
                 fontfamily = plot.symbol$fontfamily, 
#                 fontface = lattice:::chooseFace(plot.symbol$fontface, plot.symbol$font), 
                 fontface = if (is.null(plot.symbol$fontface)) plot.symbol$font else plot.symbol$fontface,
                 fontsize = trellis.par.get("fontsize")$points/2)
    plot.symbol <- trellis.par.set("plot.symbol", plot.symbol.restore)
  }
}
