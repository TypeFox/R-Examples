#' @rdname dlunif
#' @export
rlunif <- function(n, min, max, base=exp(1)) {
  if(mode(n) != "numeric")
    stop("'n' must be a non-empty numeric vector")
  if(any(missing(min), missing(max)))
    stop("'min' and 'max' not provided, without default.\n")
    ifelse(base==exp(1), return(exp(runif(n, log(min, base), log(max, base)))), return(base^(runif(n, log(min, base), log(max, base)))))
}
