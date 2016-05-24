#' 31 contrast sums from a 32 run 2^(5-0) factorial experiment on penicillin production.
#'
#' Values are arranged in decreasing order of absolute magnitude.
#' Name of the contrast effect is given as the row name of each value.
#' Daniel(1959) uses the data to illustrate the use of half-normal plots.
#' In his words: "We need, of course, some rule of inference that will help us to 
#' be objective in judging whether or not the largest effects are real."
#' 
#' \code{qqtest(penicillin[1],dist="half-normal")} will effect Daniel's plot.
#'
#'
#' @format A data frame with 31 rows and 1 variate:
#' \describe{
#'   \item{value}{value of contrast for that row}
#' }
#' @source 
#' "Use of Half-Normal Plots in Interpreting Factorial Two-Level Experiments", 
#' Cuthbert Daniel, Technometrics, Vol. 1, No. 4 (Nov., 1959), pp. 311-341.  
#' Daniel cites: Davies, O. L., Editor: Design and Analysis of Industrial Experiments, 
#'              Second Edition, Oliver and Boyd, London, and Hafner, New York, 1956 
#'              as the original. 
"penicillin"

penicillin <- data.frame(value=c(
							224, 190, 153,  93,  
							77,  64,  58,  58,
							54,  53,  53,  47,  
							39,  34,  33,  31,  
							30,  29,  28,  22,  
							21,  18,  16,  14,  
							12,   9,   7,   6,   
							4,   2,   0),
						  row.names=c(
						  "E", "A", "C", "CE",
						  "ABCDE (blocks)", "AB", "ABCD", "ACE",
						  "AD", "AC", "BC", "ACDE",
						  "BCE", "ABD", "ACD", "ABCE",
						  "DE",  "BE", "BDE", "ABE",
						  "ADE", "BCD", "BCDE", "ABDE",
						  "CDE", "D", "BD", "B",
						  "CD", "AE", "ABC"))
