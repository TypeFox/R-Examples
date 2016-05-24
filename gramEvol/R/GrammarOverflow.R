NewGrammarOverFlow <- function(...) {
  ret = list(overflow = TRUE, ...)
  class(ret) = "GrammarOverflow"
  return (ret)
}

print.GrammarOverflow <- function(x, ...) {
  print("Grammar Overflow", ...)
}

is.GrammarOverflow <- function(object) {
  return (class(object) == "GrammarOverflow")
}
