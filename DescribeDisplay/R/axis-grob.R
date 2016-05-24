#' Axes grob
#' Construct grob for axes.
#'
#' @param axes describe display object
#' @param gp arguments to the axes
#' @author Hadley Wickham \email{h.wickham@@gmail.com}
#' @keywords internal
#' @export
axesGrob <- function(axes, gp = gpar(col = "black")) {
  if (is.null(axes)) return()

  if (!is.null(axes$y)) {
    # 2d tour
    # bigaxes <- subset(as.data.frame(axes), r > 0.3)

    bigaxes <- axes[axes[, "r"] > 0.3, ]

    gTree(children = gList(
      circleGrob(
        0, 0, 1,
        default.units = "native",
        gp = gpar(fill = "transparent", col = "black")
      ),
      segmentsGrob(
        0, 0,
        axes$x, axes$y,
        default.units = "native"
      ),
      textGrob(
        bigaxes$label,
        1.1 * cos(bigaxes$theta), 1.1 * sin(bigaxes$theta),
        default.units = "native"
      )
    ), name = "axis", vp = vpPath("axes"), gp = gp)

  } else {
    # 1d tour
    n <- nrow(axes)

    gTree(children = gList(
      rectGrob(),
      linesGrob(x = unit(c(0, 0), "native"), y = unit(c(0, 1), "npc")),
      segmentsGrob(
        -1, 1:n, 1, 1:n,
        default.units = "native",
        gp = gpar(lty = 3)
      ),
      segmentsGrob(
        0, 1:n, axes$x, 1:n,
        default.units = "native",
        gp = gpar(lwd = 2)
      ),
      textGrob(
        -1:1, -1:1, -0.3,
        default.units = "native",
        just = c("centre", "top"),
        gp = gpar(cex = 0.9)
      ),
      textGrob(
        axes$label,
        1.1, 1:n,
        default.units = "native", just = c("left", "centre")
      )
    ), name = "axis", vp = vpPath("axes"), gp = gp)
  }
}

#' Axes viewport
#' Construct viewport for axes.
#'
#' @param axes describe display object
#' @param axislocation location of axes (as x and y position in npc coordinates, ie. between 0 and 1)
#' @author Hadley Wickham \email{h.wickham@@gmail.com}
#' @keywords internal
#' @export
axesViewport <- function(axes, axislocation) {
  if (is.null(axes)) return()

  if (!is.null(axes$y)) {
    # 2d tour
    viewport(
      xscale = c(-1, 1),
      yscale = c(-1, 1),
      name = "axes",
      width = 0.2, height = 0.2,
      x = axislocation[1], y = axislocation[2],
      default.units = "snpc"
    )
  } else {
    # 1d tour
    n <- nrow(axes)
    viewport(
      xscale = c(-1, 1),
      yscale = c(0, n + 1),
      name = "axes",
      width = 0.1, height = unit(n + 1, "lines"),
      x = axislocation[1], y = axislocation[2]
    )
  }
}
