##' Order a data.frame by missingness and all columns
##'
##' Completely order all rows in a data.frame.
##' 
##' @param observed a data.frame holding ordered factors in every column
##' @return the sorted order of the rows
##' @examples
##' df <- as.data.frame(matrix(c(sample.int(2, 30, replace=TRUE)), 10, 3))
##' mask <- matrix(c(sample.int(3, 30, replace=TRUE)), 10, 3) == 1
##' df[mask] <- NA
##' df[orderCompletely(df),]
orderCompletely <- function(observed) {
	observedNames <- colnames(observed)
	nacount <- sapply(observedNames, function(x) { sum(is.na(observed[,x])) })
	observedNames <- observedNames[order(nacount, decreasing=TRUE)]
	othervectorsNA <- lapply(observedNames, function(x) {!is.na(observed[,x]) })
	othervectors <- lapply(observedNames, function(x) {observed[,x] })
	args <- c(othervectorsNA, othervectors, 'na.last'=FALSE)
	do.call('order', args)
}

##' Tabulate data.frame rows
##'
##' Like \code{tabulate} but entire rows are the unit of tabulation.
##' The data.frame is not sorted, but must be sorted already.
##'
##' @param observed a sorted data.frame holding ordered factors in every column
##' @seealso \code{\link{orderCompletely}}
##' @examples
##' df <- as.data.frame(matrix(c(sample.int(2, 30, replace=TRUE)), 10, 3))
##' df <- df[orderCompletely(df),]
##' tabulateRows(df)
tabulateRows <- function(observed) {
	selectMissing <- rep(0L, nrow(observed))
	selectDefvars <- rep(0L, nrow(observed))
	threeVectors <- .Call(findIdenticalRowsData, observed,
			      selectMissing, selectDefvars, TRUE, TRUE)
	dups <- threeVectors[[1]]
	result <- rep(NA, sum(dups==1L))
	dx <- 1L
	rx <- 1L
	while (dx <= length(dups)) {
		result[rx] <- dups[dx]
		rx <- rx + 1L
		dx <- dx + dups[dx]
	}
	result
}

#' Expand summary table of patterns and frequencies
#'
#' Expand a summary table of unique response patterns to a full sized
#' data-set.
#'
#' @param tabdata An object of class \code{data.frame} with the unique response patterns and the number of frequencies
#' @param freqName Column name containing the frequencies
#' @return Returns a data frame with all the response patterns
#' @author Based on code by Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @examples
#' data(LSAT7)
#' expandDataFrame(LSAT7, freqName="freq")
expandDataFrame <- function(tabdata, freqName=NULL) {
	if (is.null(colnames(tabdata))) stop("colnames are required")

    if (missing(freqName)) {
        freqCol <- ncol(tabdata)
	warning(paste("Assuming column", colnames(tabdata)[freqCol], "contains frequencies"))
    } else {
        freqCol <- which(freqName == colnames(tabdata))
        if (length(freqCol) != 1) {
            stop(paste("Frequency column", freqName, "not found"))
        }
    }

    rows <- sum(tabdata[,freqCol])
    indexVector <- rep(NA, rows)
    rx <- 1L
    ix <- 1L
    while (rx <= nrow(tabdata)) {
        indexVector[ix:(ix + tabdata[rx,freqCol] - 1)] <- rx
	ix <- ix + tabdata[rx,freqCol]
        rx <- rx + 1L
    }
    tabdata[indexVector,-freqCol]
}

#' Compress a data frame into unique rows and frequencies
#'
#' Compress a data frame into unique rows and frequency counts.
#'
#' @param tabdata An object of class \code{data.frame}
#' @param freqColName Column name to contain the frequencies
#' @return Returns a compressed data frame
#' @examples
#' df <- as.data.frame(matrix(c(sample.int(2, 30, replace=TRUE)), 10, 3))
#' compressDataFrame(df)
compressDataFrame <- function(tabdata, freqColName="freq") {
	if (!is.na(match(freqColName, colnames(tabdata)))) {
		# Might be nice to recompress instead of stopping.
		# There might be rows to collapse due to removal
		# of columns.
		stop(paste("Frequency column", freqColName, "already appears as a column:",
			   paste(colnames(tabdata), collapse=", ")))
	}
	tabdata <- tabdata[orderCompletely(tabdata),]
	# freqs must be in numeric format, OpenMx expects integers to be factors
	freq <- as.numeric(tabulateRows(tabdata))
	tabdata <- unique(tabdata)
	tabdata <- cbind(tabdata, freq)
	colnames(tabdata)[ncol(tabdata)] <- freqColName
	tabdata
}
