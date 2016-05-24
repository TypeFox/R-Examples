
#' Create an object that describes a sankey plot
#'
#' @details
#' The node and edges data frames may contain columns that specify
#' how the plot is created. All parameters have reasonable default
#' values.
#'
#' Current list of graphical parameters for nodes:
#' \itemize{
#'   \item \code{col} Node color.
#'   \item \code{size} Node size.
#'   \item \code{x} Horizontal coordinates of the center of the node.
#'   \item \code{y} Vertical coordinates of the center of the node.
#'   \item \code{shape} Shape of the node. Possible values:
#'     \code{rectangle}, \code{point}, \code{invisible}.
#'   \item \code{lty} Lite type, see \code{par}.
#'   \item \code{srt} How to rotate the label, see \code{par}.
#'   \item \code{textcol} Label color.
#'   \item \code{label} Label text. Defaults to node name.
#'   \item \code{adjx} Horizontal adjustment of the label. See
#'     \code{adj} in the \code{par} manual.
#'   \item \code{adjy} Vertical adjustment of the label. See
#'     \code{adj} in the \code{par} manual.
#'   \item \code{boxw} Width of the node boxes.
#'   \item \code{cex} Label size multiplication factor.
#'   \item \code{top} Vertical coordinate of the top of the node.
#'   \item \code{center} Vertical coordinate of the center of the node.
#'   \item \code{bottom} Vertical coordinate of the bottom of the node.
#'   \item \code{pos} Position of the text label, see \code{par}.
#'   \item \code{textx} Horizontal position of the text label.
#'   \item \code{texty} Vertical position of the text label.
#' }
#'
#' Current list of graphical parameters for edges:
#' \itemize{
#'   \item \code{colorstyle} Whether the to use a solid color (\code{col}),
#'     or \code{gradient} to plot the edges. The color of a gradient
#'     edges is between the colors of the nodes.
#'   \item \code{curvestyle} Edge style, \code{sin} for sinusoid curves,
#'     \code{line} for straight lines.
#'   \item \code{col} Edge color, for edges with solid colors.
#'   \item \code{weight} Edge weight. Determines the width of the edges.
#' }
#'
#' @param nodes A data frame of nodes on the plot, and possibly
#'   their visual style. The first column must be the ids of the
#'   nodes. If this argument is \code{NULL}, then the ids of the
#'   nodes are determined from \code{edges}.
#' @param edges A data frame of the edges. The first two columns
#'   must be node ids, and they define the edges. The rest of the columns
#'   contain the visual style of the edges.
#' @param y How to calculate vertical coordinates of nodes, if they
#'   are not given in the input. \code{optimal} tries to minimize edge
#'   crossings, \code{simple} simply packs nodes in the order they are
#'   given, from bottom to top.
#' @param break_edges Whether to plot each edge as two segments,
#'   or a single one. Sometimes two segment plots look better.
#' @param gravity Whether to push the nodes to the top, to the bottom
#'   or to the center, within a column.
#' @return A \code{sankey} object that can be plotted via the
#'   \code{\link{sankey}} function.x
#'
#' @importFrom simplegraph graph
#' @export

make_sankey <- function(
  nodes = NULL, edges, y = c("optimal", "simple"), break_edges = FALSE,
  gravity = c("center", "top", "bottom")) {

  y <- match.arg(y)
  gravity <- match.arg(gravity)
  
  if (is.null(nodes)) {
    nodes <- data.frame(
      stringsAsFactors = FALSE,
      id = sort(unique(c(edges[,1], edges[,2])))
    )
  }

  ## Easy ones first, breaking edges depend on these, so we need
  ## to add them early

  nodes[["col"]]     <- nodes[["col"]]     %||% color_nodes(nodes, edges)
  nodes[["shape"]]   <- nodes[["shape"]]   %||% "rectangle"
  nodes[["lty"]]     <- nodes[["lty"]]     %||% 1
  nodes[["srt"]]     <- nodes[["srt"]]     %||% 0
  nodes[["textcol"]] <- nodes[["textcol"]] %||% "black"
  nodes[["label"]]   <- nodes[["label"]]   %||% nodes[,1]
  nodes[["adjx"]]    <- nodes[["adjx"]]    %||% 1/2
  nodes[["adjy"]]    <- nodes[["adjy"]]    %||% 1
  nodes[["boxw"]]    <- nodes[["boxw"]]    %||% 0.2
  nodes[["cex"]]     <- nodes[["cex"]]     %||% 0.7

  edges[["colorstyle"]] <- edges[["colorstyle"]] %||% "gradient"
  edges[["curvestyle"]] <- edges[["curvestyle"]] %||% "sin"
  edges[["col"]]        <- edges[["col"]]        %||% color_edges(nodes, edges)
  edges[["weight"]]     <- edges[["weight"]]     %||% 1

  nodes[["size"]]    <- nodes[["size"]]    %||% optimize_sizes(nodes, edges)
  nodes[["x"]]       <- nodes[["x"]]       %||% optimize_x(nodes, edges)

  ## We can break the edges now
  if (break_edges) {
    ne <- do_break_edges(nodes, edges)
    nodes <- ne$nodes
    edges <- ne$edges
  }

  if (null_or_any_na(nodes[["y"]])      ||
      null_or_any_na(nodes[["top"]])    ||
      null_or_any_na(nodes[["center"]]) ||
      null_or_any_na(nodes[["bottom"]])) {
    nodes <- optimize_y(nodes, edges, mode = y, gravity = gravity)
  }

  if (null_or_any_na(nodes[["pos"]])    ||
      null_or_any_na(nodes[["textx"]])  ||
      null_or_any_na(nodes[["texty"]])) {
    nodes <- optimize_pos(nodes, edges)
  }

  ## Reorder edges in the order of node ids, so that edges
  ## coming from the same node do not cross
  node_ids <- nodes[,1]
  node_order <- base::order(nodes$x, nodes$center)
  edges <- edges[ base::order(node_order[match(edges[,1], node_ids)],
                              node_order[match(edges[,2], node_ids)]), ]

  res <- graph(nodes, edges)
  class(res) <- c("sankey", class(res))
  res
}

color_edges <- function(nodes, edges) {
  "#99d8c9"
}

color_nodes <- function(nodes, edges) {
  "#2ca25f"
}

#' @importFrom simplegraph predecessors successors

optimize_sizes <- function(nodes, edges) {

  sgraph <- graph(nodes, edges)

  lefts  <- vapply(predecessors(sgraph), length, 1L)
  rights <- vapply(successors(sgraph), length, 1L)

  pmax(pmax(lefts, rights), 1)
}

optimize_pos <- function(nodes, edges) {

  ## Middle nodes are at (center, below)
  ## First column nodes are at (right, center)
  ## Last column is (left, center)

  first <- nodes$x == min(nodes$x)
  last  <- nodes$x == max(nodes$x)

  if (is.null(nodes$pos)) {
    nodes$pos <- rep(1, nrow(nodes))
    nodes$pos[first] <- 2
    nodes$pos[last ] <- 4
  }

  if (is.null(nodes$textx) || is.null(nodes$texty)) {
    nodes$textx <- nodes$x
    nodes$texty <- nodes$top

    nodes$textx [first] <- nodes$x[first] - nodes$boxw[first] / 2
    nodes$texty [first] <- nodes$center[first]

    nodes$textx [last]  <- nodes$x[last] + nodes$boxw[last] / 2
    nodes$texty [last]  <- nodes$center[last]
  }

  nodes
}
