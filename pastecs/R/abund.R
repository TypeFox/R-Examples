"abund" <-
function(x, f = 0.2) {
	call <- match.call()
	# Rem: we could decide to store the initial data into res$data
	# To free memory, we will just store a call to these data
	# The drawback is that initial data should not be modified
	# between 'abund' and 'extract.abund'!!!
	Data <- deparse(substitute(x))
	x <- as.data.frame(x)	# We want to be sure to work on a data frame!
	q <- nrow(x)
	# Percentage of non-null values in each series
	percnonull <- apply((x != 0), 2, sum)/q*100
	percnull <- 100 - percnonull
	# Log of total number of individuals counted in each series + 1, in percent of the most abundant one
	nbrind <- log(apply(x, 2, sum))
	# If there were some columns full of zeros, log() return -Inf => reset them to 0
	nbrind[is.infinite(nbrind)] <- 0
	percind <- nbrind/max(nbrind)*100
	# We adjust the value of f if necessary
	if (!is.numeric(f)) 
		stop("abund: f must be numeric")
	if (f <= 0) {
		if (f != 0) warning("f must not be smaller than zero. Was adjusted to 0.00001")
		f <- 0.00001
	}
	if (f >= 1) {
		if (f != 1) warning("f must not be larger than one. Was adjusted to 0.99999")
		f <- 0.99999
	}
	# Perform the sorting of descriptors
	srt <- f*percind + (1 - f)*percnull		# Calculate criterion for sort of 'local abundance'
	srtlist <- sort.list(srt)				# Sort in ascending order
	# Sort descriptors accordingly and rescale them for cumsum calculation
	srt <- srt[srtlist]
	percindsrt <- percind[srtlist]
	pi <- (percindsrt-min(percindsrt))/(max(percindsrt)-min(percindsrt))*100
	percnonullsrt <- percnonull[srtlist]
	pn <- (percnonullsrt-min(percnonullsrt))/(max(percnonullsrt)-min(percnonullsrt))*100
	cumsm <- cumsum(pi-pn)
	# Rescale cumsm between 0 and 100
	cumsm <- ((cumsm-min(cumsm))/(max(cumsm)-min(cumsm)))*100
	names(cumsm) <- names(percindsrt)
	names(srtlist) <- names(percindsrt)
	names(srt) <- names(percindsrt)
	res <- list(data=Data, vr=srtlist, sort=srt, cumsum=cumsm, p.log.ind=percindsrt, p.nonull=percnonullsrt, f=f, call=call)		# Create a list containing the result
	class(res) <- "abund"						# and turn it into an 'abund' object
	res											# Return the result
}
