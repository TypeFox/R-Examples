Bolt <- function(FUNC, listen=0, boltID=0){
	x <- list()
	if(!is.function(FUNC)){
		stop("Please provide a processing function")
	}
	x$name <- as.character(substitute(FUNC))
	x$func <- FUNC
	x$listen <- listen
	x$id <- boltID
	class(x) <- "Bolt"
	x
}

# Check, internal function
is.Bolt <- function(x){
	ifelse( class(x) == "Bolt", TRUE, FALSE)
}
