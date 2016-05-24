kNanosecond <- 1
kMicrosecond <- kNanosecond * 1000
kMillisecond <- kMicrosecond * 1000
kSecond <- kMillisecond * 1000
kMinute <- kSecond * 60
kHour <- kMinute * 60

formatDurationSmall <- function(ns) {
    sizes <- c('ns', 'us', 'ms', 's')
    e <- ifelse(ns == 0, NA, floor(log(ns, 1000)))
    suffix <- ifelse(ns == 0, '', sizes[e+1])
    prefix <- ifelse(ns == 0, '0', sprintf("%g", ns/(1000 ^ floor(e))))
    paste(prefix, suffix, sep="")
}

formatDurationHour <- function(ns) {
	sprintf("%.0fh%s", ns/kHour, formatDurationMinute(ns %% kHour))
}

formatDurationMinute <- function(ns) {
	ifelse(ns > kHour,
		formatDurationHour(ns),
		sprintf("%.0fm%ss", ns/kMinute,
			format((ns %% kMinute) / kSecond, digits=9, scientific=F)))
}

formatDuration <- function(ns) {
	prefix <- ifelse(ns < 0, "-", "")
	ns <- abs(ns)
	paste(prefix, ifelse(ns > kMinute,
			formatDurationMinute(ns),
			formatDurationSmall(ns)), sep="")
}

formatNanoseconds <- formatDuration

formatMicroseconds <- function(us) {
	formatNanoseconds(us * kMicrosecond)
}

formatMilliseconds <- function(ms) {
	formatNanoseconds(ms * kMillisecond)
}

formatSeconds <- function(s) {
	formatNanoseconds(s * kSecond)
}
