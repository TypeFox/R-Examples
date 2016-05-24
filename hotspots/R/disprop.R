disprop <-
function(z) {
if (!inherits(z, "hotspots"))
	stop("use only with \"hotspots\" objects")
dp <- list(positive = (z$x-median(z$data))/
	(z$positive.cut - median(z$data)),
negative = (z$x-median(z$data))/
	(z$negative.cut - median(z$data)))
if (z$tail == "positive") dp$negative <- NULL
if (z$tail == "negative") dp$positive <- NULL
dp}

