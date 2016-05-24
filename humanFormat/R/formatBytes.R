formatAnyBytes <- function(b, base, fmt, sizes) {
	e <- ifelse(b == 0, NA, floor(log(b, base)))
	suffix <- ifelse(b == 0, "", sizes[e + 1])
	prefix <- ifelse(b < base, b, sprintf(fmt, b/(base^floor(e))))
	paste(prefix, suffix, sep="")
}

formatIECBytes <- function(b, fmt="%.2f") {
	formatAnyBytes(b, 1024, fmt, c("", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB", "ZiB", "YiB"))
}

formatSIBytes <- function(b, fmt="%.2f") {
	formatAnyBytes(b, 1000, fmt, c("", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"))
}

formatBytes <- formatSIBytes
