"partial.corr" <-
function(vars, cond) {
  y <- data.matrix(vars)
  cor(residuals(lsfit(y=y, x=data.matrix(cond),
                      yname=dimnames(y)[[2]])))
}
