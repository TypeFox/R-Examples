"alpha.proxy" <-
function (weight=.2, vol.man=.2, vol.bench=.2, vol.other=.2, cor.man=.2,
	cor.bench=.2, plot.it=TRUE, transpose=FALSE, ...)
{
	fun.copyright <- "Placed in the public domain 2003-2012 by Burns Statistics Ltd."
        fun.version <- "alpha.proxy 002"


	# check ranges for possible bad scaling

	if(any(weight <=0 | weight >= 1)) stop("bad value(s) for weight")
	if(any(cor.man < -1 | cor.man > 1)) stop("bad value(s) for cor.man")
	if(any(cor.bench < -1 | cor.bench > 1)) 
		stop("bad value(s) for cor.bench")
	if(any(vol.man <= 0)) stop("all vol.man values must be positive")
	if(any(vol.bench <= 0)) stop("all vol.bench values must be positive")
	if(any(vol.other <= 0)) stop("all vol.other values must be positive")
	if(all(vol.man >= 1)) {
		warning(paste("large vol.man values, did you mistakenly",
			"give values in percent?"))
	}
	if(all(vol.bench >= 1))  {
		warning(paste("large vol.bench values, did you mistakenly",
			"give values in percent?"))
	}
	if(all(vol.other >= 1))  {
		warning(paste("large vol.other values, did you mistakenly",
			"give values in percent?"))
	}

	# get organized

	sizes <- numeric(6)
	names(sizes) <- c("weight", "vol.man", "vol.bench", "vol.other",
		"cor.man", "cor.bench")
	sizes["weight"] <- length(weight)
	sizes["vol.man"] <- length(vol.man)
	sizes["vol.bench"] <- length(vol.bench)
	sizes["vol.other"] <- length(vol.other)
	sizes["cor.man"] <- length(cor.man)
	sizes["cor.bench"] <- length(cor.bench)

	if(any(sizes == 0)) stop("zero length input(s)")
	if(sum(sizes > 1) > 2) stop("more than two inputs longer than 1")
	twovecs <- sum(sizes > 1) == 2

	# do computation

	if(!twovecs) {
		ans <- -weight * (vol.man^2 - vol.bench^2) - 2 *
			(1 - weight) * vol.other * (vol.man * cor.man -
			vol.bench * cor.bench)
		return(10000 * ans)
	}

	count <- 1
	longsizes <- sizes[sizes != 1]
	z <- array(NA, longsizes)

	for(i in names(longsizes)) assign(i, sort(get(i)))

	# nested for loops sacrifice efficiency for convenience
	for(i.cb in 1:sizes["cor.bench"]) {
	for(i.cm in 1:sizes["cor.man"]) {
	for(i.vo in 1:sizes["vol.other"]) {
	for(i.vb in 1:sizes["vol.bench"]) {
	for(i.vm in 1:sizes["vol.man"]) {
	for(i.w in 1:sizes["weight"]) {
		z[count] <- -weight[i.w] * (vol.man[i.vm]^2 - 
			vol.bench[i.vb]^2) - 2 * (1 - weight[i.w]) * 
			vol.other[i.vo] * (vol.man[i.vm] * cor.man[i.cm] -
			vol.bench[i.vb] * cor.bench[i.cb])
		count <- count + 1
	}}}}}}
	z <- 10000 * z

	# actually do something

	if(transpose) {
		z <- t(z)
		longsizes <- rev(longsizes)
	}
	ans <- list(x = eval(as.name(names(longsizes[1]))),
		y = eval(as.name(names(longsizes[2]))), z = z,
		call = deparse(match.call()))
	if(plot.it && twovecs) {
		the.labs <- c(weight="Weight", vol.man="Manager Volatility",
			vol.bench="Benchmark Volatility",
			vol.other="Volatility of the Rest",
			cor.man="Correlation of Manager with the Rest",
			cor.bench="Correlation of Benchmark with the Rest")
		filled.contour(ans, xlab=the.labs[names(longsizes)[1]],
			ylab=the.labs[names(longsizes)[2]], 
			plot.axis={axis(1); axis(2); contour(ans, add=TRUE)},
			...)
		invisible(ans)
	} else {
		ans
	}
}

