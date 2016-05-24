
optimal_edge_order <- function(nodes, edges) {

  sel <- function(x, attr) nodes[ match(x, nodes[,1]), attr]

  left_x  <- sel(edges[,1], "x")
  right_x <- sel(edges[,2], "x")

  left_y  <- sel(edges[,1], "center")
  right_y <- sel(edges[,2], "center")

  base::order(-left_x, right_y, right_x)
}

#' @importFrom simplegraph vertices edges strength

draw.edges <- function(x, nsteps = 50) {
  # for each node, we need to to store the position of the current slot, on
  # the right and on the left

  nodes <- vertices(x)
  edges <- edges(x)

  nodes$lpos <- nodes$center - strength(x, mode = "in")  / 2
  nodes$rpos <- nodes$center - strength(x, mode = "out") / 2

  edge_order <- optimal_edge_order(nodes, edges)

  for (i in edge_order) {

    n1 <- edges[i, 1]
    n2 <- edges[i, 2]

    sel <- function(node, attr) nodes[ nodes[,1] == node, attr]

    curveseg(
      sel(n1, "x") + sel(n1, "boxw") / 2,
      sel(n2, "x") - sel(n2, "boxw") / 2,
      sel(n1, "rpos"),
      sel(n2, "lpos"),
      colorstyle = edges$colorstyle[i],
      grad = c(sel(n1, "col"), sel(n2, "col")),
      width = edges$weight[i],
      col = edges$col[i],
      nsteps = nsteps,
      curvestyle = edges$curvestyle[i]
    )

    nodes[nodes[,1] == n1, "rpos"] <- sel(n1, "rpos") + edges$weight[i]
    nodes[nodes[,1] == n2, "lpos"] <- sel(n2, "lpos") + edges$weight[i]
  }
}

#' @importFrom graphics points rect text

draw.nodes <- function(x, width = 0.2) {

  nodes <- vertices(x)

  for (n in seq_len(nrow(nodes))) {

    if (nodes$shape[n] == "invisible") next

    if (nodes$shape[n] == "point") {
      points(nodes$x[n], nodes$center[n], pch = 19, col = nodes$col[n])

    } else if (nodes$shape[n] == "rectangle") {
      rect(
        nodes$x[n] - nodes$boxw[n] / 2, nodes$bottom[n],
        nodes$x[n] + nodes$boxw[n] / 2, nodes$top[n],
        lty = nodes$lty[n], col = nodes$col[n]
      )
    }

    text(
      nodes$textx[n],
      nodes$texty[n],
      nodes$label[n],
      col = nodes$textcol[n],
      srt = nodes$srt[n],
      pos = nodes$pos[n],
      adj = c(nodes$adjx[n], nodes$adjy[n]),
      cex = nodes$cex[n],
      offset = 0.2,
      xpd = NA
    )
  }
}

#' @rdname sankey
#' @method plot sankey
#' @export

plot.sankey <- function(x, ...) sankey(x, ...)

#' Draw a sankey plot
#'
#' @param x The plot, created via \code{\link{make_sankey}}.
#' @param mar Margin of the plot, see \code{mar} in the \code{par}
#'   manual.
#' @param ... Additional arguments, ignored currently.
#' @return Nothing.
#'
#' @export
#' @importFrom graphics par plot.new
#' @importFrom grDevices dev.hold dev.flush

sankey <- function(x, mar = c(0, 5, 0, 5) + 0.2, ...) {

  plot.new()

  V <- vertices(x)
  E <- edges(x)

  xrange <- range(V$x, na.rm = TRUE)
  xlim <- xrange + (xrange[2] - xrange[1]) * c(-0.1, 0.1)
  yrange <- range(V[, c("bottom", "top")], na.rm = TRUE)
  ylim <- yrange + (yrange[2] - yrange[1]) * c(-0.1, 0.1)

  par(mar = mar)
  par(usr = c(xlim, ylim))

  dev.hold()
  on.exit(dev.flush())

  draw.edges(x, nsteps = 50)
  draw.nodes(x)

  invisible()
}
