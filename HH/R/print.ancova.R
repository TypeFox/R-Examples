"print.ancova" <-
function(x, ...) {
  print(anova(x, ...))
  print(attr(x,"trellis"))
  invisible(x)
}

