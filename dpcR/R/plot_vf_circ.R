plot_vf_circ <- function(x, y, radius, bg) {
  symbols(x = x, y = rep(y, length(x)), circles = rep(radius, length(x)), 
          bg = bg, fg = "black", inches = FALSE, add = TRUE)
}