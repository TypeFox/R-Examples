origin <- structure(0, class = c("POSIXt", "POSIXct"))

POSIXt.units <- function(x=NULL, ...) {
	availables.units <- c('second', 'minute', 'hour', 'day',
			      'week', 'month', 'year', 'AD')
	if(is.null( x ))
		return(factor(availables.units, availables.units,
			      ordered = TRUE) ) else
		return (factor (x, availables.units, ordered=TRUE) )
	      # return (factor (availables.units[c(1:4, 6:7)],
	      # availables.units, ordered = TRUE) ) else
}

setMethod ('Ops', c('numeric', 'ANY'),
	function (e1, e2) do.call(.Generic, list(e1=e1, e2=as.numeric(e2))))
setMethod ('Ops', c('ANY', 'numeric'),
	function (e1, e2) do.call(.Generic, list(e1=as.numeric(e1), e2=e2)))
