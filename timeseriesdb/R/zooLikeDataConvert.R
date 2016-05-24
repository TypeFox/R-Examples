#' Zoo like Date Conversion
#' 
#' This function is taken from the zoo package. It is basically the
#' S3 method as.Date.numeric of the package zoo. It is used to turn
#' 2005.75 (3rd quarter of 2005) like date formats into dates
#' like 2005-07-01.
#'
#' @param x object of class ts
#' @param offset numeric defaults to 0. See the zoo package for more information.
#' @param ... optional arguments. 
#' @rdname zooLikeDateConvert
#' @name zooLikeDateconvert
#' 
#' @author Achim Zeileis, Gabor Grothendieck, Jeffrey A. Ryan,
#' Felix Andrews
#' @export
zooLikeDateConvert <- function (x, offset = 0, ...) 
{
  time.x <- unclass(time(x)) + offset
  if (frequency(x) == 1) 
    as.Date(paste(time.x, 1, 1, sep = "-"))
  else if (frequency(x) == 4) 
    as.Date(paste((time.x + 0.001)%/%1, 3 * (cycle(x) - 1) + 
                    1, 1, sep = "-"))
  else if (frequency(x) == 12) 
    as.Date(paste((time.x + 0.001)%/%1, cycle(x), 1, sep = "-"))
  else stop("unable to convert ts time to Date class")
}
