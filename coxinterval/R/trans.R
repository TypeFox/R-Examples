trans <- function(from, to)
{
  x <- cbind(table(from, to), sort(unique(from)), sort(unique(to)))
  state1 <- x[apply(x[, 1:2] == 0, 1, any), 3]
  state0 <- x[x[, 3] != state1, 3]
  state2 <- x[x[, 4] != state1, 4]
  x <- cbind(from, to)
  attr(x, "states") <- c(state0, state1, state2)
  cl <- match.call()
  types <- cl[c(1, match(c("from", "to"), names(cl), nomatch = 0))]
  types[[1]] <- as.name("strata")
  attr(x, "types") <- levels(eval(types, parent.frame()))
  x
}
