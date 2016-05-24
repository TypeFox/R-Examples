# create the constraint that w_i >= w_j
ordinalConstraint <- function(n, i, j) {
  a <- rep(0, n)
  a[i] <- -1
  a[j] <- 1
  list(constr=t(a), rhs=c(0), dir=c("<="))
}

# create the constraint that w_i <= x
upperBoundConstraint <- function(n, i, x) {
  a <- rep(0, n)
  a[i] <- 1
  list(constr=t(a), rhs=c(x), dir=c("<="))
}

# create the constraint that w_i >= x
lowerBoundConstraint <- function(n, i, x) {
  a <- rep(0, n)
  a[i] <- -1
  list(constr=t(a), rhs=c(-x), dir=c("<="))
}

# create the constraint that w_i/w_j <= x
upperRatioConstraint <- function(n, i, j, x) {
  a <- rep(0, n)
  a[i] <- 1
  a[j] <- -x
  list(constr=t(a), rhs=c(0), dir=c("<="))
}

# create the constraint that x <= w_i/w_j
lowerRatioConstraint <- function(n, i, j, x) {
  a <- rep(0, n)
  a[i] <- -1
  a[j] <- x
  list(constr=t(a), rhs=c(0), dir=c("<="))
}

# create the constraint that x = w_i/w_j
exactRatioConstraint <- function(n, i, j, x) {
  a <- rep(0, n)
  a[i] <- -1
  a[j] <- x
  list(constr=t(a), rhs=c(0), dir=c("="))
}

# create the n-simplex
simplexConstraints <- function(n) {
  list(constr = rbind(rep(1, n), -diag(n)),
       dir = c('=', rep('<=', n)),
       rhs = c(1, rep(0, n)))
}

# merge a list of constraints
mergeConstraints <- function(...) {
  lst <- list(...)
  if (length(lst) == 1 && is.list(lst[[1]])) {
    lst <- lst[[1]]
  }
  list(
    constr=do.call("rbind", lapply(lst, function(c) { c$constr })),
    rhs=unlist(lapply(lst, function(c) { c$rhs })),
    dir=unlist(lapply(lst, function(c) { c$dir }))
  )
}

# filter a set of constraints
filterConstraints <- function(constr, sel) {
  list(constr = constr[['constr']][sel, , drop=FALSE],
       rhs = constr[['rhs']][sel],
       dir = constr[['dir']][sel])
}
