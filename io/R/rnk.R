#' @method qread rnk
#' @export
qread.rnk <- function(file, type, rm.forced=FALSE, ...) {
	d <- read.table(file, sep="\t", row.names=NULL, header=FALSE,
		colClasses = c("character", "numeric"), ...);
	x <- d[,2];
	names(x) <- d[,1];

	x
}

#' @method qwrite rnk
#' @export
qwrite.rnk <- function(x, file, type, value.var=NULL, ...) {
	if (is.matrix(x)) {
		stopifnot(!is.null(rownames(x)));
		if (is.null(value.var)) {
			stopifnot(ncol(x) == 1);
			value.var <- 1;
		}
		d <- data.frame(
			feature = rownames(x),
			value = x[,1]
		);
	} else if (is.data.frame(x)) {
		if (is.null(value.var)) {
			# x should be already pre-formatted
			stopifnot(ncol(x) == 2);
			d <- x;
		} else {
			stopifnot(!is.null(rownames(x)));
			d <- data.frame(
				feature = rownames(x),
				value = x[, value.var]
			);
		}
	} else if (is.numeric(x)) {
		stopifnot(!is.null(names(x)));
		d <- data.frame(
			feature = names(x),
			value = x
		);
	} else {
		stop("Unsupported data type");
	}

	stopifnot(is.numeric(d[,2]));
	stopifnot(all(!is.na(d[,2])));

	write.table(d, file, row.names=FALSE, col.names=FALSE, quote=FALSE,
		sep="\t", ...)
}
