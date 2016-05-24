#' @method [<- lfactor
#' @export
"[<-.lfactor" <- function(x,i, value) {
	cc <- value %in% llevels(x)
	ll <- llevels(x)
	lapply(1:length(llevels(x)), function(j) {
		value[value==ll[j]] <- levels(x)[j]
		value <<- value
  })
	NextMethod()
}
