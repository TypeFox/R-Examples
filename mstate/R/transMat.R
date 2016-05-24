transMat <- function(x, names) {
  ## transMat:  produce transition matrix for use in package 'mstate'.
  ## Arguments:
  ## x:  List of possible transitions.
  ##     x[[i]] consists of a vector of state numbers
  ##     reachable from state i.
  ## names: Character vector of state names, having the same length
  ##        as x.
  ## Example:  States 1 and 2 are reachable from one
  ##     another.  State 3 is absorbing, reachable from
  ##     states 1 and 2.
  ##     transMat( x = list( c(2, 3), c(1, 3), c() ) )
  if ( !is.list(x) ) stop("x must be a list")
  ns <- length(x) ## number of states
  tmat <- matrix(NA, nrow = ns, ncol = ns) ## transition matrix
  if ( missing(names) ) {
    if ( !is.null( base::names(x) ) ) {
      namesList <- list(from = base::names(x), to = base::names(x))
    } else {
      namesList <- list(from = paste("State", seq(nrow(tmat))),
                        to = paste("State", seq(nrow(tmat))))
    }
  } else {
    if ( length(names) != ns ) stop("length of 'names' must equal length of 'x'")
    namesList <- list(from = names, to = names)
  }
  idxmat <- cbind(unlist(lapply(seq(ns),
                                function(i, y){
                                  rep(i, length(y[[i]]))}, x)),
                  unlist(x))
  if ( max(idxmat) > ns )
    stop("Largest state in transition list exceeds number of states")
  tmat[idxmat] <- seq(nrow(idxmat))
  dimnames(tmat) <- namesList
  tmat
}
