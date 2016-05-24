#' @include precintcon.read.data.r
NULL

#' @name read.data
#' @aliases precintcon.read.data read.data 
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' 
#' @title Load a precipitation series 
#' @description Load a file with a daily or monthly precipitation series. 
#' @usage read.data(file, sep = ",", dec = ".", header = TRUE, na.value = NA) 
#' @param file a string containing the file path.
#' @param sep the character applied for delimited columns. (Default value: ",")
#' @param dec the character applied for defined decimal point. (Default value: ".")
#' @param header a logical value defining whether the first line of the file refers to column names. (Default value: TRUE)
#' @param na.value the value used for representing missing values. (Default value: NA)
#' @return 
#'   A \code{data frame} containing a representation of the
#'   data in the \code{file}. The \code{file} is addressed as precintcon.daily or 
#'   precintcon.monthly depending of its structure.
#'   
#'   The file should contains three columns when loading monthly series 
#'   and thirty three columns when loading daily series.
#'   
#'   The first columns refers to years and the second one refers to months.
#'   When dealing with daily datasets, the thirty one remaining columns refers 
#'   to the amount of precipitation in the days of the months.
#'   Otherwise, the remaining column refers to the amount of precipitation in each month.   
#' @seealso 
#'   \code{\link{daily}}
#'   \code{\link{monthly}}
#'   \code{\link{read.table}}
#'   \code{\link{read.csv}}
#'   \code{\link{read.csv2}}
#' @examples 
#' ##
#' # Loading a serie on Windows
#' \dontrun{d1 <- read.data("C:\PRECINTCON\203040.csv", sep = ";", dec = ".", header = TRUE)}
#' 
#' ##
#' # Loading a serie on Unix-like
#' \dontrun{d1 <- read.data("/home/precintcon/203040.csv", sep = ";", dec = ".", header = TRUE)}
#' @keywords import read data read table file precipitation 
#' @export
read.data <- function(file, sep = ",", dec = ".", header = TRUE, na.value = NA)
   precintcon.read.data(file = file, sep = sep, dec = dec, header = header, na.value = na.value)