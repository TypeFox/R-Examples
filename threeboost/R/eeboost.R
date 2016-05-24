#' EEBoost
#' 
#' Alias for ThrEEBoost (which defaults to a threshold value of 1).
#' @export
#' @param ... Arguments to \code{\link{threeboost}}.
#' @seealso \code{\link{threeboost}}
eeboost <- function(...) { threeboost(...) }