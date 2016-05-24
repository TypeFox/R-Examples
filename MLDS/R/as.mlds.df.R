`df2mlds.df` <- function(d, st) {
	.Deprecated("as.mlds.df", "MLDS")
	if (missing(st)) st <- sort(unique(unlist(d[, -1])))
	names(d) <- c("resp", paste("S", 1:4, sep = ""))
	attr(d, "stimulus") <- st
	invord <- d$S1 > d$S4
	attr(d, "invord") <- invord
	class(d) <- c("mlds.df", "data.frame")
	d
}

`as.mlds.df` <- function(d, ...)
	UseMethod("as.mlds.df")
	
#`as.mlds.df.default` <- function(d, ...)
#    NextMethod("as.mlds.df", d, ...)

`as.mlds.df.data.frame` <- function(d, st, ...) {
	cl <- if (any(d == -2)) 3 else 4
	cls <- if(cl == 4) "mlds.df" else "mlbs.df"
	if (ncol(d) > 5) {
		stim1 <- ifelse(rowSums(d[, -1]) == 0, 0, 1)
		ix.mat <- cbind(stim.1 = stim1, d[, -1])
		dd <- t(apply(ix.mat, 1, function(x) which(x != 0)))
		dd <- as.data.frame(cbind(d[, 1], dd))
		names(dd) <- c("resp", paste("S", 1:cl, sep = ""))
		d <- dd
		if (missing(st)) st <- sort(unique(unlist(d[, -1])))
		attr(d, "stimulus") <- st
		d
		} else {
	if (missing(st)) st <- sort(unique(unlist(d[, -1])))
	names(d) <- c("resp", paste("S", 1:cl, sep = ""))
	attr(d, "stimulus") <- st
	invord <- d[[2]] > d[[length(d)]]
	attr(d, "invord") <- invord
	class(d) <- c(cls, "data.frame")
	d
	}
}

`as.mlbs.df` <- function(d, ...)
	UseMethod("as.mlbs.df")
	
#`as.mlbs.df.default` <- function(d, ...)
#    NextMethod("as.mlbs.df", d, ...)
	
`as.mlbs.df.data.frame` <- function(d, st, ...) {
	if (missing(st))	
		st <- sort(unique(unlist(d[, -1])))
	attr(d, "stimulus") <- st
	names(d) <- c("resp", paste("S", 1:3, sep = ""))
	invord <- d$S1 > d$S3
    attr(d, "invord") <- invord
    class(d) <- c("mlbs.df", "data.frame")
    d
}      	