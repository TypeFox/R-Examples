
optimize_y <- function(nodes, edges, mode = c("optimal", "simple"),
                       gravity) {

  mode <- match.arg(mode)

  if (mode == "simple") {
    optimize_y_simple(nodes, edges, gravity = gravity)

  } else if (mode == "optimal") {
    optimize_y_optim(nodes, edges, gravity = gravity)
  }
}


optimize_y_simple <- function(nodes, edges,
                              gravity = c("top", "bottom", "center")) {

  gravity <- match.arg(gravity)

  ## 10 percent of max total node size at a level
  interstop <- 0.3 * max(tapply(nodes$size, nodes$x, sum))

  nodes$center <- nodes$top <- nodes$bottom <- NA_real_

  if (gravity == "top") {
    optimize_y_simple_top(nodes, interstop)
  } else if (gravity == "bottom") {
    optimize_y_simple_bottom(nodes, interstop)
  } else {
    optimize_y_simple_center(nodes, interstop)
  }
}


optimize_y_simple_top <- function(nodes, interstop) {

  for (pos in sort(unique(nodes$x))) {
    cur_y <- 0
    nodes_here <- rev(which(nodes$x == pos))

    for (node in nodes_here) {

      if (! is.null(nodes$y) && ! is.na(nodes$y[node])) {
        nodes$center[node] <- nodes$y[node]
        nodes$top[node]    <- nodes$y[node] - nodes$size[node] / 2
        nodes$bottom[node] <- nodes$y[node] + nodes$size[node] / 2

      } else {
        nodes$bottom[node] <- cur_y
        nodes$center[node] <- cur_y - nodes$size[node] / 2
        nodes$top   [node] <- cur_y - nodes$size[node]
        cur_y <- cur_y - nodes$size[node] - interstop
      }
    }
  }

  nodes
}


optimize_y_simple_bottom <- function(nodes, interstop) {

  for (pos in sort(unique(nodes$x))) {
    cur_y <- 0
    nodes_here <- which(nodes$x == pos)

    for (node in nodes_here) {

      if (! is.null(nodes$y) && ! is.na(nodes$y[node])) {
        nodes$center[node] <- nodes$y[node]
        nodes$top[node]    <- nodes$y[node] - nodes$size[node] / 2
        nodes$bottom[node] <- nodes$y[node] + nodes$size[node] / 2

      } else {
        nodes$top   [node] <- cur_y
        nodes$center[node] <- cur_y + nodes$size[node] / 2
        nodes$bottom[node] <- cur_y + nodes$size[node]
        cur_y <- cur_y + nodes$size[node] + interstop
      }
    }
  }

  nodes
}


optimize_y_simple_center <- function(nodes, interstop) {

  for (pos in sort(unique(nodes$x))) {
    cur_y <- 0
    nodes_here <- which(nodes$x == pos)

    for (node in nodes_here) {

      if (! is.null(nodes$y) && ! is.na(nodes$y[node])) {
        nodes$center[node] <- nodes$y[node]
        nodes$top[node]    <- nodes$y[node] - nodes$size[node] / 2
        nodes$bottom[node] <- nodes$y[node] + nodes$size[node] / 2

      } else {
        nodes$bottom[node] <- cur_y
        nodes$center[node] <- cur_y - nodes$size[node] / 2
        nodes$top   [node] <- cur_y - nodes$size[node]
        cur_y <- cur_y - nodes$size[node] - interstop
      }
    }
  }

  ylim <- range(nodes$bottom, nodes$top)
  dy <- ylim[2] - ylim[1]

  for (pos in sort(unique(nodes$x))) {
    nodes_here <- which(nodes$x == pos)
    ylim_here <- range(nodes[nodes_here, c("bottom", "top")])
    dy_here   <- ylim_here[2] - ylim_here[1]

    ## Only center if it is not centered
    if ( sum(ylim_here) != sum(ylim) ) {
      nodes[nodes_here, c("bottom", "center", "top")] <-
        nodes[nodes_here, c("bottom", "center", "top")] - (dy - dy_here) / 2
    }
  }

  nodes
}
