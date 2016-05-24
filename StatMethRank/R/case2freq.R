#' Convert raw ranking data(case form) to ranking data with rank frequencies
#'
#' Convert raw ranking data(case form) to ranking data with rank frequencies
#'
#' @param x A data frame or a matrix(case form), each row is a rank.
#' @return A data frame, each row contains a rank and the corresponding frequency.
#' @export
#' @author Li Qinglong <liqinglong0830@@163.com>
#' @examples
#' data(APA)
#' cases = freq2case(APA, freq.col = 1)
#' freqs = case2freq(cases)
case2freq <- function (x) 
{
	nCol = dim(x)[2]
	DF = as.data.frame(table(x), stringAsFactors = TRUE)
	DF = DF[DF[, nCol + 1] != 0, ]
	DF = as.data.frame(DF)
	row.names(DF) = 1:nrow(DF)
	return(DF)
}

