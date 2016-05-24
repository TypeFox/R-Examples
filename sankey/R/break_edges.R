
## We can assume that nodes and edges have color at
## this point. We add labels here, if not present,
## to hide the extra nodes.

do_break_edges <- function(nodes, edges) {

  sel <- function(id, attr) nodes[match(id, nodes[,1]), attr]

  to_break <- which(sel(edges[,2], "x") - sel(edges[,1], "x") > 1)

  if (length(to_break) == 0) return(list(nodes = nodes, edges = edges))

  ## Color of new nodes is the mean of the two incident nodes
  col1 <- sel(edges[to_break, 1], "col")
  col2 <- sel(edges[to_break, 2], "col")
  col <- mean_colors(col1, col2)

  new_nodes <- data.frame(
    stringsAsFactors = FALSE,
    id = make.unique(
      paste(edges[to_break, 1], sep = "-", edges[to_break, 2]), sep = "_"
    ),
    x = (sel(edges[to_break, 1], "x") + sel(edges[to_break, 2], "x")) / 2,
    label = "",
    size = 1,
    shape = "invisible",
    boxw = 0,
    col = col
  )
  names(new_nodes)[1] <- names(nodes)[1]

  edges2 <- edges[to_break, ]
  edges[to_break, 2] <- new_nodes[,1]
  edges2[,1] <- new_nodes[,1]

  edges <- rbind(edges, edges2)
  nodes <- merge(nodes, new_nodes, all = TRUE)

  list(nodes = nodes, edges = edges)
}

mean_colors <- function(col1, col2) {
  vapply(seq_along(col1), FUN.VALUE = "", function(i) {
    mrgb <- rowMeans(cbind(col2rgb(col1[i]), col2rgb(col2[i])))
    do.call(rgb, as.list(mrgb / 255))
  })
}
