latex.wprof <-
function(y, label="", caption="", scale=c(1, 1), ...) {
  Z <- getzeta(y)
  latex(Z, label, caption, scale)
}
