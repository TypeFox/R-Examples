getResponse <- function(formula) {
  tt <- terms(formula)
  vars <- as.character(attr(tt, "variables"))[-1] 
  response <- attr(tt, "response")
  vars[response] 
}