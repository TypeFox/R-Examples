
nodes_edges <- function() {

  nodes <- data.frame(
    stringsAsFactors = FALSE,
    name = c("add_mid_points_fun", "autoscale", "calcpos",
      "calcsizes2", "draw.edges", "draw.nodes", "orderConnections",
      "plot.sankey", "sankey", "color_ramp_palette_alpha", "ypos_present",
      "curveseg", "points", "rect", "text", "plot.new", "checkedges",
      "par", "dev.hold", "dev.flush", "strwidth")
  )

  edges <- data.frame(
    stringsAsFactors = FALSE,
    from = c("add_mid_points_fun", "autoscale", "draw.edges",
      "draw.nodes", "draw.nodes", "draw.nodes", "plot.sankey", "sankey",
      "sankey", "sankey", "sankey", "sankey", "sankey", "sankey", "sankey",
      "sankey", "sankey", "sankey", "sankey", "sankey", "sankey", "sankey"),
    to = c("color_ramp_palette_alpha", "ypos_present", "curveseg",
      "points", "rect", "text", "sankey", "plot.new", "checkedges",
      "orderConnections", "ypos_present", "add_mid_points_fun", "autoscale",
      "calcsizes2", "calcpos", "par", "par", "dev.hold", "dev.flush",
      "strwidth", "draw.edges", "draw.nodes")
  )

  list(nodes = nodes, edges = edges)
}
