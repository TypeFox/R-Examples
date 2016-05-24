#' @export
precintcon.pn.analysis <- function(
   object, 
   interval = 30, 
   scale    = "a"
) {
	
	if (!is.element(scale, c("w", "m", "s", "a", "d")))
		stop("Invalid scale parameter. It should be either \"w\" (weak), \"m\" (month), \"s\" (seasonal), \"a\" (annual), or \"d\" (decade)")

	if (interval < 1)
		stop("invalid interval. It should be greater than 0")
		
	l <- precintcon.pn(object, interval=interval, scale=scale)		
	
	return(l)		
}
