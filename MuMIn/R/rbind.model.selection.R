`rbind.model.selection` <-
function (..., deparse.level = 1, make.row.names = TRUE) {
	allargs <- list(...)
	n <- length(allargs)
	if(n == 1L) return(allargs[[1L]])

	if(!all(vapply(allargs, inherits, FALSE, "model.selection")))
		stop("need all \"model.selection\" objects")

	### XXX: This modifies original objects!!!
	allargs <- lapply(allargs, "class<-", "data.frame")
	## ... reverting to original (?) class on exit:
	on.exit(lapply(allargs, "class<-", c("model.selection", "data.frame")))

	allitemsidentical <- function(x) all(vapply(x[-1L], identical, FALSE, x[[1L]]))

	if(!allitemsidentical(lapply(lapply(allargs, attr, "rank"), attr, "call")))
		stop("tables are not ranked by the same IC")
	if(!allitemsidentical(lapply(allargs, "attr", "nobs")))
		stop("models are fitted to different number of observations")

	.combine <-
	function(x, y, pos, len = length(y)) {
		if(is.factor(x) || is.factor(y)) {
			if(is.factor(x)) {
				if(!is.factor(y)) y <- factor(y)
			} else if(is.factor(y)) x <- factor(x)
			alllev <- unique(c(levels(x), levels(y)))
			x <- factor(x, levels = alllev, labels = alllev)
		}
		x[pos:(pos + len - 1L)] <- y
		x
	}

	ct <- unname(lapply(allargs, attr, "column.types"))
	vct <- unlist(ct, recursive = FALSE)
	vct <- vct[order(as.integer(unlist(ct)))]

	#vct <- vct[order(as.integer(unlist(ct)), unlist(lapply(ct, seq_along)))]
	vct <- vct[!duplicated(names(vct))]
	# TODO: check mismatch in column.types
	nm <- names(vct)

	rval <- as.data.frame(array(NA, dim = c(sum(sapply(allargs, nrow)), length(nm)),
								dimnames = list(NULL, nm)))
	row1 <- 1L
	for(z in allargs) {
		n <- nrow(z)
		nmz <- nm[nm %in% names(z)]
		for(j in nmz) rval[, j] <- .combine(rval[, j], z[, j], row1, n)
		row1 <- row1 + n
	}

	newattr <- list(column.types = vct)
	for(i in c("model.calls", "coefTables"))
		newattr[[i]] <- unlist(lapply(allargs, attr, i), recursive = FALSE, use.names = FALSE)
	k <- c("rank", "nobs")
	newattr[k] <- attributes(allargs[[1L]])[k]

	tmp <- lapply(allargs, attr, "terms")
	newattr[["terms"]] <- structure(unique(unlist(tmp, recursive = FALSE, use.names = FALSE)),
			  interceptLabel = unique(unlist(lapply(tmp, attr, "interceptLabel"))))


	for(i in names(newattr)) attr(rval, i) <- newattr[[i]]
	class(rval) <- c("model.selection", "data.frame")
	if(make.row.names) {
		rn1 <- rep(names(allargs), sapply(allargs, nrow))
		rn1[i] <- paste0(rn1[i <- rn1 != ""], ".")
		rlabs <- paste0(rn1, unlist(lapply(allargs, rownames)))
		if(anyDuplicated(rlabs))
			rlabs <- make.unique(as.character(rlabs), sep = "")
	} else {
		rlabs <- as.character(1L:nrow(rval))
	}
	rownames(rval) <- rlabs

	o <- order(rval[, names(vct)[vct == "ic"]])
	rval <- rval[o, recalc.delta = TRUE]
	attr(rval, "merged-order") <- o
	rval
}

`merge.model.selection` <-
function (x, y, suffixes = c(".x", ".y"), ...) {
	rval <- rbind(x, y, make.row.names = FALSE)
	if (!is.null(suffixes)) row.names(rval) <-
		c(paste0(row.names(x), suffixes[1L]),
            paste0(row.names(y), suffixes[2L]))[attr(rval, "merged-order")]
	attr(rval, "merged-order") <- NULL
	rval
}
