#' Convert frequency table to raw data 
#'
#' Convert frequency table to raw data(case form). This is a useful function to 
#' convert the data written by Michael Friendly.
#'
#' @param x A table object, or a data frame in frequency form containing ranks and one
#' numeric variable representing the frequency for that rank.
#' @param var.names A list of variable names for the factors, if you wish to override those already in
#' the table.
#' @param freq.name The name of the frequency variable in the table.
#' @param freq.col The column index of the frequencies.
#' @param ... Other arguments passed down to type.convert.
#' @return A data frame containing the factors in the table and as many observations as are represented by the
#' total of the freq variable.
#' @export
#' @author Michael Friendly, edited by Li Qinglong <liqinglong0830@@163.com>
#' @examples
#' data(APA)
#' cases = freq2case(APA, freq.col = 1)
#' freqs = case2freq(cases)
#' @references Posted on R-Help, Jan 20, 2009.


freq2case <- function (x, var.names = NULL, freq.name = "Freq", freq.col = NULL, ...) 
{
	# Convert frequency table to raw data 
	# source code : expand.dft function of vcdExtra package, Author: Michael Friendly
	# convert frequency table to raw data 
    # Editted by Li Qinglong
    
    if (inherits(x, "table")) 
        x <- as.data.frame.table(x, responseName = freq.name)
    if (is.null(freq.col))
        freq.col <- which(colnames(x) == freq.name)
    if (length(freq.col) == 0) 
        stop(paste(sQuote("freq"), "not found in column names"))
    DF <- sapply(1:nrow(x), function(i) x[rep(i, each = x[i, 
        freq.col]), ], simplify = FALSE)
    DF <- do.call("rbind", DF)[, -freq.col]
    for (i in 1:ncol(DF)) {
        DF[[i]] <- type.convert(as.character(DF[[i]]), ...)
    }
    rownames(DF) <- NULL
    if (!is.null(var.names)) {
        if (length(var.names) < dim(DF)[2]) {
            stop(paste("Too few", sQuote("var.names"), "given."))
        }
        else if (length(var.names) > dim(DF)[2]) {
            stop(paste("Too many", sQuote("var.names"), "given."))
        }
        else {
            names(DF) <- var.names
        }
    }
    return(DF)
}