#' @rdname dlunif
#' @export
plunif <- function(q, min, max, base=exp(1))  {
  if(mode(q) != "numeric")
    stop("'q' must be a non-empty numeric vector")
  if(any(missing(min), missing(max)))
    stop("'min' and 'max' not provided, without default.\n")
  return(1/(log(max,base)-log(min,base)) * (log(q,base)-log(min,base)))
}
