#' Get quantile contour
#'
#' @param x Output from \code{\link[MASS]{kde2d}}.
#' @param alpha The quantile level.
get_quantile_contour <- function(x, alpha = 0.8) {
  zdens <- rev(sort(x$z))
  Czdens <- cumsum(zdens)
  Czdens <- (Czdens/Czdens[length(zdens)])
  crit.val <- zdens[max(which(Czdens<=alpha))]
  b.full=contourLines(x,levels=crit.val)
  list(x = b.full[[1]]$x, y = b.full[[1]]$y)
}

#' Custom bandwidth
#'
#' Based on \code{bandwidth.nrd} from \pkg{MASS}. This version takes the
#' absolute value of \code{var} to avoid errors.
#'
#' @param x A numeric vector
custom_bw <- function(x) {
    r <- quantile(x, c(0.25, 0.75))
    h <- (r[2] - r[1])/1.34
    4 * 1.06 * min(sqrt(abs(var(x))), h) * length(x)^(-1/5)
}

#' Add a kernel density polygon
#'
#' @param x x values
#' @param y y values
#' @param col Colour to add polygon with. Will be made into two levels of
#'   opacity.
#' @param lwd lwd Line width
#' @param alpha A numeric vector of length 2 that gives the confidence levels
#'   for the two kernel density polygons.
#' @param add_poly Add polygons?
#' @param add_pts Logical: should points be added?
add_dens_polygon <- function(x, y, col, lwd = 1.7, alpha = c(0.25, 0.75), add_pts = FALSE, add_poly = TRUE) {
  if(add_poly) {
    x_bw <- custom_bw(na.omit(x))
    y_bw <- custom_bw(na.omit(y))
    k <- get_quantile_contour(MASS::kde2d(x,y, h = c(x_bw, y_bw)), alpha = 0.75)
    polygon(k$x, k$y, border = col, col = paste(col, "40", sep = ""), lwd = lwd)
    k <- get_quantile_contour(MASS::kde2d(x,y, h = c(x_bw, y_bw)), alpha = 0.25)
    polygon(k$x, k$y, border = col, col = paste(col, "80", sep = ""), lwd = lwd)
  }
  if(add_pts) points(x, y, pch = 21, col = paste(col, "60", sep = ""), cex = 0.7)
}

#' Get the efficient frontier from mean and variance values
#'
#' @param m A vector of mean values
#' @param v A vector of variance values
get_efficient_frontier <- function(m, v) {
  d <- data.frame(m = m, v = v)
  ef_front <- d[chull(d$v, d$m), ]

  ef_front <- ef_front[order(ef_front$m), ]
  ef_front_top <- ef_front[which.min(ef_front$v):nrow(ef_front), ]
  ef_front_bottom <- ef_front[1:(which.min(ef_front$v) - 1), ]

  # order by risk
  ef_front_top <- ef_front_top[order(ef_front_top$v), ]
  ef_front_bottom <- ef_front_bottom[order(-ef_front_bottom$v), ]

  # remove where return is coming back down or back up:
  ef_front_top <- ef_front_top[1:which.max(ef_front_top$m),]
  ef_front_bottom <- ef_front_bottom[which.min(ef_front_bottom$m):nrow(ef_front_bottom),]
  return(rbind(ef_front_bottom, ef_front_top))
}

#' Plot conservation plans in mean-variance space
#'
#' This makes a mean-variance plot of the portfolio output. It can take care of:
#' plotting the individual portfolios, adding 2D kernel density polygons at two
#' quantile levels, and adding an efficient frontier.
#'
#' @param plans_mv The \code{plans_mv} element of the output from
#'   \code{\link{run_cons_plans}}.
#' @param plans_name A character vector of what to label each conservation
#'   plan.
#' @param cols Colours for the conservation plan polygons.
#' @param xlim X limits
#' @param ylim Y limits
#' @param add_pts Logical: add the points?
#' @param add_all_efs Logical: add efficient frontiers?
#' @param x_axis Logical: add x axis?
#' @param y_axis Logical: add y axis?
#' @param add_legend Logical: add y legend?
#' @param legend_pos A character string to pass to
#' \code{\link[graphics]{legend}} denoting the position of the legend.
#' @param w_show If \code{"all"} then all plans will be shown. If a numeric
#'   vector, then those plans will be shown. E.g. \code{c(1, 3)} will only show
#'   the first and third plans.
#' @param xlab X axis label.
#' @param ylab Y axis label.
#' @param ... Anything else to pass to \code{\link[graphics]{plot.default}}.
#' @param add_poly Add the kernal smoother quantile polygons?
#' @return A plot. Also, the x and y limits are returned invisibly as a list.
#'   This makes it easy to make the first plot and then save those x and y
#'   limits to fix them in subsequent (multipanel) plots.
#' @export
plot_cons_plans <- function(plans_mv, plans_name, cols, xlim = NULL,
  ylim = NULL, add_pts = TRUE, add_all_efs = FALSE, x_axis = TRUE,
  y_axis = TRUE, add_legend = TRUE, legend_pos = "topright",
  w_show = "all", xlab = "Variance", ylab =
  "Mean", add_poly = TRUE, ...) {

  if(w_show[1] == "all") w_show <- seq_along(plans_name)

  plans_mv <- plyr::llply(plans_mv, function(x) {
    x_out <- na.omit(x)
    if(nrow(x_out) < nrow(x)) warning("Some simulations are not plotted due to NAs.")
    x_out
  })

  if(is.null(xlim)) {
    lims <- plyr::ldply(plans_mv, function(x) data.frame(x.max =
        max(x$v, na.rm = TRUE), x.min = min(x$v, na.rm = TRUE), y.max
        = max(x$m, na.rm = TRUE), y.min = min(x$m, na.rm = TRUE)))
    xlim = c(min(lims$x.min, na.rm = TRUE), max(lims$x.max, na.rm = TRUE)*0.8)
    ylim = c(min(lims$y.min, na.rm = TRUE), max(lims$y.max, na.rm = TRUE))
  }

  plot.default(1, 1, xlim = xlim, ylim = ylim, type = "n", xlab = xlab, ylab = ylab,
    axes = FALSE, ...)
    box(col = "grey50")
    if(x_axis) axis(1, col=  "grey50", at = pretty(axTicks(1), n = 4))
    if(y_axis) axis(2, col=  "grey50", at = pretty(axTicks(2), n = 4))

    for(i in w_show) {
      add_dens_polygon(plans_mv[[i]]$v, plans_mv[[i]]$m, col = cols[i],
        add_pts = add_pts, add_poly = add_poly)
    }
    if(add_legend) {
    legend(legend_pos, legend = plans_name[w_show], fill =
      #paste(cols[w_show], "95", sep = ""), bty = "n", cex = 1.0)
      paste(cols[w_show], sep = ""), bty = "n", cex = 1.0)
    }

    mv_all <- do.call("rbind", plans_mv)

    ef_front_pts <- with(mv_all, get_efficient_frontier(m = m, v = v))
    with(ef_front_pts, lines(v, m, col = "grey50", lwd = 2.2))

    if(add_all_efs) {
      for(i in w_show) {
        ef_front_pts_i <- with(plans_mv[[i]], get_efficient_frontier(m =
            m, v = v))
        with(ef_front_pts_i, lines(v, m, col = cols[i], lwd = 2.2))
      }
    }
    invisible(list(xlim = xlim, ylim = ylim))
}
