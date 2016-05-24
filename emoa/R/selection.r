##
## selection.r - Selection schemes for EAs
##
## All selection methods should have the same signature:
##
##   selection(values, n, ...)
##

##' Selection strategies for EMOA.
##'
##' The currently implemented strategies are nondominated sorting
##' followed by either hypervolume contribution or crowding distance
##' based ranking. Both of these implementations are currently
##' limited to selecting a single individual for replacement.
##'
##' @param values Matrix of function values.
##' @param n      Number of individuals to select for replacement.
##' @param ...    Optional parameters passed to
##'   \code{\link{hypervolume_contribution}}. 
##' 
##' @title Selection strategies
##' @aliases nds_hv_selection nds_cd_selection
##'
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
##' @keywords optimize nonlinear
##' @export
nds_hv_selection <- function(values, n=1, ...) {  
  #stopifnot(n == 1)
  sel <- which(is_maximally_dominated(values))
  
  ## Identify individual which gets replaced:
  if (length(sel) == 1) {
    sel
  } else {
    contrib <- if (length(sel) == ncol(values)) {
        hypervolume_contribution(values, ...)
    } else {
       hypervolume_contribution(values[,sel], ...)
     }
    sel[which.min(contrib)]
  }
}

##' @export
##' @rdname nds_hv_selection
nds_cd_selection <- function(values, n=1, ...) {
  #stopifnot(n == 1)
  N <- ncol(values)
  k <- N - n
  ranks <- nds_rank(values)
  sel <- rep(FALSE, N)
  cr <- 0
  while(sum(sel) < k) {
    cr <- cr + 1
    sel[ranks == cr] <- TRUE
  }

  if (sum(sel) != k) {
    nelim <- sum(sel) - k
    dist <- crowding_distance(values[,ranks == cr])
    cdr <- rank(dist, ties.method="random")
    s <- which(ranks == cr)[cdr <= nelim]
    sel[s] <- FALSE
  }
  which(sel == FALSE)
}
