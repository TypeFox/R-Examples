#' Draw dd plot
#' Draw a complete describe display.
#'
#' If you want to layout multiple dd plots on the same page, you can
#' use \code{\link[grid]{grid.layout}}.  If you need even more control,
#' set \code{draw = FALSE} and then \code{\link[grid]{grid.draw}} the
#' resulting grob yourself.
#'
#' This function reads a number of options directly out of the
#' descripedisplay datastructure.  See the examples for ways to use
#' these.
#'
#' @param x dd object to plot
#' @param ... (unused)
#' @param draw draw plot, or just return grob
#' @param axislocation location of axes (as x and y position in npc coordinates, ie. between 0 and 1)
#' @param size size of plot as a proportion of the total display area (set to 1 for printed out)
#' @param axisgp color of the axis
#' @param background.color color of in the background of the plot
#' @return frame grob containing all panels, note that this does not contain the title or border
#' @examples
#' plot(dd_example("xyplot"))
#' plot(dd_example("tour1d"))
#' plot(dd_example("tour2d"))
#' @author Hadley Wickham \email{h.wickham@@gmail.com}
#' @keywords internal
#' @method plot dd
#' @export
plot.dd <- function(
  x,
  ...,
  draw = TRUE,
  axislocation = c(0.1, 0.1),
  size = 0.9,
  axisgp = gpar(col = "black"),
  background.color = "grey90"
) {
  d <- x$dim
  layout <- grid.layout(nrow = d[1], ncol = d[2])
  panels <- frameGrob(layout = layout)

  for (i in 1:x$nplot) {
    panels <- placeGrob(panels,
      ddpanelGrob(
        x$plots[[i]],
        axislocation = axislocation,
        axis.gp = axisgp,
        background.color = background.color
      ),
      col = (i - 1) %/% d[1] + 1, row = (i - 1) %% d[1] + 1
    )
  }

  if (!is.null(x$title) && nchar(x$title) != 0) {
    pg <- frameGrob(grid.layout(nrow = 2, ncol = 1))
    pg <- packGrob(
      pg,
      textGrob(
        x$title,
        gp = gpar(cex = 1.3)
      ),
      row = 1, height = unit(2, "lines")
    )
    pg <- packGrob(pg, panels, row = 2)
  } else {
    pg <- panels
  }

  if (draw) {
    grid.newpage()
    pushViewport(viewport(width = size, height = size))
    grid.draw(pg)
  }

  invisible(panels)
}

#' Plot a dd plot
#' Convenient method to draw a single panel.
#'
#' This is mainly used for bug testing so that you can pull out a single
#' panel quickly and easily.
#'
#' @param x object to plot
#' @param ... (not used)
#' @param axislocation location of axes (as x and y position in npc coordinates, ie. between 0 and 1)
#' @param axis.gp frame grob containing all panels, note that this does not contain the title or border
#' @param background.color color of in the background of the plot
#' @author Hadley Wickham \email{h.wickham@@gmail.com}
#' @keywords hplot
#' @method plot ddplot
#' @export
#' @examples
#' scatmat <- dd_example("scatmat")
#' plot(scatmat)
#' plot(scatmat$plots[[1]])
#' plot(scatmat$plots[[3]])
#' plot(scatmat$plots[[4]])
plot.ddplot <- function(
  x,
  ...,
  axislocation = c(0.1, 0.1),
  axis.gp = gpar(col = "black"),
  background.color = "grey90"
) {
  grid.newpage()
  grob <- ddpanelGrob(
    x,
    axislocation = axislocation,
    axis.gp = axis.gp,
    background.color = background.color
  )
  grid.draw(grob)
}


#' Panel grob
#' Construct grob for single panel.
#'
#' @param panel describe display object
#' @param axis location, x and y position
#' @param axisgp frame grob containing all panels, note that this does not contain the title or border
#' @param background.color color of in the background of the plot
#' @author Hadley Wickham \email{h.wickham@@gmail.com}
#' @keywords internal
ddpanelGrob <- function(
  panel,
  axislocation = c(0.1, 0.1),
  axis.gp = gpar(col = "black"),
  background.color = "grey90"
) {
  points <- panel$points
  edges <- panel$edges

  axes <- dd_tour_axes(panel)
  axesVp <- axesViewport(axes, axislocation)
  grobs <- list(
    rectGrob(gp = gpar(col = "grey", fill = background.color))
  )

  if (!is.null(edges)) {
    grobs <- append(
      grobs,
      list(
        segmentsGrob(
          edges$src.x, edges$src.y,
          edges$dest.x, edges$dest.y,
          default.units = "native",
          gp = gpar(
            lwd = edges$lwd,
            col = as.character(edges$col)
          )
        )
      )
    )
  }

  if (is.null(panel$showPoints) || panel$showPoints) {
    grobs <- append(
      grobs,
      list(
        pointsGrob(
          points$x, points$y,
          pch = points$pch,
          gp = gpar(
            col = as.character(points$col)
          ),
          size = unit(points$cex, "char")
        )
      )
    )
  }

  if (!is.null(panel$labels)) {
    labels <- panel$labels
    grobs <- append(
      grobs,
      list(
        textGrob(
          as.character(labels$label),
          labels$x, labels$y,
          default.units = "native",
          hjust = labels$left, vjust = labels$top
        )
      )
    )
  }

  grobs <- append(grobs,  list(
    textGrob(panel$params$xlab %||% "", 0.99, 0.01,
      just = c("right", "bottom")),
    textGrob(panel$params$ylab %||% "", 0.01, 0.99, just = c("left", "top")),
    axesGrob(axes, gp = axis.gp)
  ))

  if (length(panel$params$label) == 1)
    grobs <- append(grobs, list(textGrob(panel$params$label %||% "",
     0.5, 0.01, just = c("centre", "bottom"))))

  if (!is.null(panel$drawlines) && panel$drawlines) {
    grobs <- append(
      grobs,
      list(
        segmentsGrob(
          points$x,
          panel$baseline,
          points$x, points$y,
          default.units = "native",
          gp = gpar(
            col = as.character(points$col)
          )
        )
      )
    )
  }


  gTree(
    children = do.call(gList, grobs),
    vp = dataViewport(
      xscale = panel$xscale,
      yscale = panel$yscale,
      clip = "on"),
    childrenvp = axesVp
  )
}
