round.POSIXlt <- function(x, digits = c("mins", "5mins", "10mins", "15mins", "quarter hours", "30mins", "half hours", "hours")){
	units <- digits
	# code found on https://stat.ethz.ch/pipermail/r-help/2012-June/315336.html
	if(is.numeric(units)){
		r <- 60*units[1]
	}
	else{
		units <- match.arg(digits)
		r <- switch(units,
			"mins" = 60,
			"5mins" = 60*5,
			"10mins" = 60*10,
			"15mins"=, "quarter hours" = 60*15,
			"30mins"=, "half hours" = 60*30,
			"hours" = 60*60
		)
	}
	H <- as.integer(format(x, "%H"))
	M <- as.integer(format(x, "%M"))
	S <- as.integer(format(x, "%S"))
	D <- format(x, "%Y-%m-%d")
	secs <- 3600*H + 60*M + S
	as.POSIXct(round(secs/r)*r, origin=D, tz=attributes(x)$tzone)
}