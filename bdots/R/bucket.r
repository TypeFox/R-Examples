bucket <- function(sig, times, ylim) {
	buckets <- bucket.calc(sig, times)
	bucket.plot(buckets, ylim)
	invisible(buckets)
}
