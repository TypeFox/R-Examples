switchllevels <- function(x) {
  factor(x, levels=levels(x), labels=llevels(x))
}
