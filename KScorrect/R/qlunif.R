#' @rdname dlunif
#' @export
qlunif <- function(p, min, max, base=exp(1)) {
  if(mode(p) != "numeric")
    stop("'p' must be a non-empty numeric vector")
  if(any(missing(min), missing(max)))
    stop("'min' and 'max' not provided, without default.\n")
  return(base ^ (log(min, base) + (log(max, base) - log(min, base)) * p))
  }
