#' Time series with one of each violation to be used as test and example sample.
#'
#' @docType data
#' @details
#' It was built with the following code:\cr\cr
#' code{values = cbind(values=c(17,\cr
#' 39, # rule 1, point 2 red\cr
#' 11,12,10,5,16,8,5,15,14, # rule 2, point 11 yellow3\cr
#' 21,\cr
#' 17,18,19,21,22,25, # rule 3, point 18 green\cr
#' 6,\cr
#' 7,18,9,22,18,21,16,21,18,20,18,24,18,23, # rule 4, point 33 magenta\cr
#' 31,15,31, # rule 5, point 36 blue\cr
#' 10,\cr
#' 25,24,9,24,26, # rule 6, point 42 orange\cr
#' 17,18,17,20,18,11,14,19,15,16,18,13,22,20,10, # rule 7, point 57 brown\cr
#' 7,25,24,8,24,8,25,24, #rule 8, point 65 cyan\cr
#' 22,5,21,16,12,11,33,17,15,13,22,13,11,8,23,5,10,6,10,21,5,9,11,20,8,23,14,19,5,12,17,7,15,14,9))\cr
#' dates = as.character(seq(from=as.Date("1/1/2015", "\%d/\%m/\%Y"), length.out = 100, by = "days"))\cr\cr
#' series <- data.frame(values, dates)\cr\cr
#' series <- series[sample(nrow(series)),]}
#' @format A data frame with 100 rows and 2 variables:
#' \describe{
#'   \item{values}{100 values ranging from 5 to 31}
#'   \item{dates}{100 dates from 1-jan-2015 to 10-apr-2015 to give the values a time order}
#' }
"spcTimeSeries"
