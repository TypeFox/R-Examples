
##' groups
##' 
##' @keywords internal
##' @param groups groups
##' @param n.groups n.groups
##' @param n 
groups <- function(groups, n.groups = NULL, n) {
  ng <- length(groups)
  if (is.null(n.groups) | !is.list(n.groups)) {
    n.groups <- list(n.groups)
  }
  n.groups <- rep(n.groups, ng)
  
  for (i in 1:ng) {
    lg <- length(groups[[i]])
    if (!is.null(n.groups[[i]]))
      n.groups[[i]] <- rep(n.groups[[i]], length = lg)
    if (is.null(n.groups[[i]])) {
      n.groups[[i]] <- rep(floor(n/lg), lg)
    }
    if (sum(n.groups[[i]]) != n) {
      n.groups[[i]][length(n.groups[[i]])] <- n.groups[[i]][length(n.groups[[i]])] + n - sum(n.groups[[i]])
    }
  }
  results <- list(groups, n.groups)
  results
}



##' ngroups
##' 
##' @keywords internal
##' @param group group
##' @param n.group n.group
##' @param n n
ngroups <- function(group, n.group = NULL, n) {
  ng <- length(group)
  pos.group <- c(1, 1 + cumsum(n.group))[1:ng]
  data.frame(group, pos.group, n.group, stringsAsFactors = FALSE)
}


##' linegroup
##' 
##' @keywords internal
##' @param group group
##' @param n.group n.group
linegroup <- function(group, n.group) {
  res <- list()
  for (i in 1:length(group)) {
    tmp <- NULL
    for (j in 1:length(group[[i]])) {
      tmp <- c(tmp, expand(group[[i]][j], 1, n.group[[i]][j], what = ""))
    }
    res <- c(res, list(tmp))
  }
  res
}
