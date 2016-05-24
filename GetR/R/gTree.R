gTree <- function(formula, data = list(), type = "once") {
  require(party)
  ctree(formula, data, controls = ctree_control(mincriterion = 0.1, stump = TRUE))
}