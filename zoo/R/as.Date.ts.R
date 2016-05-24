# as.Date.numeric <- function (x, origin = "1970-01-01", ...) 
#  as.Date(origin, ...) + x

as.Date.ts <- function(x, offset = 0, ...) {
   time.x <- unclass(time(x)) + offset
   if (frequency(x) == 1)
	as.Date(paste(time.x, 1, 1, sep = "-"))
   else if (frequency(x) == 4)
	as.Date(paste((time.x + .001) %/% 1, 3*(cycle(x)-1)+1, 1, sep = "-"))
   else if (frequency(x) == 12)
	as.Date(paste((time.x + .001) %/% 1, cycle(x), 1, sep = "-"))
   else
	stop("unable to convert ts time to Date class")
}

