latex.incidence <-
function(y, label="", caption="", scale=c(1, 1), ...) {
  C <- incidence2cover(y)
  latex(C, label, caption, scale)
}
