
optimize_y_optim <- function(nodes, edges, gravity) {

  ## Starting state
  nodes <- optimize_y_simple(nodes, edges, gravity = gravity)
  nodes$y <- nodes$center

  ## But reorder nodes according to their x and then y coords
  nodes <- nodes[ base::order(nodes$x, nodes$y), ]

  ## And we also rewrite the y coordinates to 1:k
  nodes <- set_integer_y(nodes)

  xpos <- sort(unique(nodes$x))
  for (pos in xpos) nodes <- bubble(nodes, edges, pos)

  ## Need to run it again, to fix vertical space between nodes
  ## But first reorder the nodes in good order
  nodes <- nodes[ base::order(nodes$x, nodes$y), ]
  nodes$y <- NULL
  optimize_y_simple(nodes, edges, gravity = gravity)
}

set_integer_y <- function(nodes) {
  xpos <- sort(unique(nodes$x))
  for (pos in xpos) {
    num <- sum(nodes$x == pos)
    nodes[nodes$x == pos, ]$y <- seq_len(num)
  }
  nodes
}

switch_nodes <- function(nodes, x, y1, y2) {
  idx1 <- which(nodes$x == x & nodes$y == y1)
  idx2 <- which(nodes$x == x & nodes$y == y2)
  nodes[idx1,]$y <- y2
  nodes[idx2,]$y <- y1
  nodes
}


## Take the nodes at a certain position and put them
## into optimal order. This is a double loop, we bring
## up each node from the bottom, individually. It is similar
## to bubble sort.
##
## We identify the nodes with their y positions.
## This must be unique.
bubble <- function(nodes, edges, pos) {
  ypos <- nodes$y[nodes$x == pos]
  for (pos2 in rev(ypos)) nodes <- bring_up(nodes, edges, pos, pos2)
  for (pos2 in ypos) nodes <- bring_down(nodes, edges, pos, pos2)
  nodes
}

crossing_edge <- function(nodes, edges, e1, e2) {
  t1 <- nodes[ edges[e1, 1] == nodes[,1], ]$y
  h1 <- nodes[ edges[e1, 2] == nodes[,1], ]$y
  t2 <- nodes[ edges[e2, 1] == nodes[,1], ]$y
  h2 <- nodes[ edges[e2, 2] == nodes[,1], ]$y
  sign(t1 - t2) + sign(h1 - h2) == 0
}

crossing_edges <- function(nodes, edges, eset1, eset2) {
  crossing <- 0
  for (e1 in eset1) {
    for (e2 in eset2) {
      crossing <- crossing + crossing_edge(nodes, edges, e1, e2)
    }
  }
  crossing
}

## Count the crossing edges for incoming and outgoing
## edges of nodes at (x,y1) and (x,y2)
eval_node_pair <- function(nodes, edges, x, y1, y2) {
  idx1 <- which(nodes$x == x & nodes$y == y1)
  idx2 <- which(nodes$x == x & nodes$y == y2)
  edges_in1 <- which(edges[,2] == nodes[idx1, 1])
  edges_in2 <- which(edges[,2] == nodes[idx2, 1])
  edges_ou1 <- which(edges[,1] == nodes[idx1, 1])
  edges_ou2 <- which(edges[,1] == nodes[idx2, 1])

  crossing_edges(nodes, edges, edges_in1, edges_in2) +
    crossing_edges(nodes, edges, edges_ou1, edges_ou2)
}

## Take node at (x,y) and bring it up as much as possible.
bring_up <- function(nodes, edges, x, y) {
  while (y > 1) {
    nodes <- try_switch_nodes(nodes, edges, x, y, y - 1)
    y <- y - 1
  }
  nodes
}

bring_down <- function(nodes, edges, x, y) {
  maxy <- max(nodes$y[nodes$x == x])
  while (y < maxy) {
    nodes <- try_switch_nodes(nodes, edges, x, y, y + 1)
    y <- y + 1
  }
  nodes
}

try_switch_nodes <- function(nodes, edges, x, y1, y2) {

  befor <- eval_node_pair(nodes, edges, x, y1, y2)
  nodes <- switch_nodes(nodes, x, y1, y2)
  after <- eval_node_pair(nodes, edges, x, y1, y2)
  if (befor < after) {
    nodes <- switch_nodes(nodes, x, y1, y2)
  }
  nodes
}
