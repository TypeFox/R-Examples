##' A method of the generic function by for objects of class \'trackdata\'
##' 
##' A given function 'FUN' is applied to the data corresponding to each segment
##' of data.
##' 
##' trapply() applies a function iteratively to each segment of a trackdata
##' object without the need for using a for-loop. It can be used to calculate,
##' for example, the mean value of the data values of each segment separately.
##' Any function that can be applied sensibly to trackdata[j]\$data where j is
##' a segment number can be used as the fun argument to trapply(). It is also
##' possible to write your own function and use trapply() to apply it
##' separately to each segment. Care needs to be taken in using trapply() in
##' the following two ways. Firstly, the argument simplify=T should only be set
##' if it can be guaranteed that a vector of the same length or matrix of the
##' same number of rows as the number of segments in the trackdata object is
##' returned. For example, simplify=T can be used in calculating the mean per
##' segment of a trackdata object, because there will only be one value (the
##' mean) per segment. However, simplify should be set to F in calculating the
##' range because here two values are returned per segment. Similarly use
##' simplify=F n smoothing the data in which the number of values returned per
##' segment is different.  Secondly, trapply() only applies a function to a
##' single parameter; the function can be used to apply to a function to
##' multi-parameter trackdata such as F1-F4, but then the function needs to be
##' put inside apply() - see examples below.
##' 
##' @param trackdata a track data object
##' @param fun a function that is applied to each segment
##' @param \dots arguments of the function fun
##' @param simplify simplify = TRUE , output is a matrix; simplify = FALSE a
##' list is returned
##' @param returntrack returntrack = FALSE , return a trackdata object
##' @return list or vector or matrix
##' @author Jonathan Harrington
##' @seealso \code{\link{apply}}
##' @keywords methods
##' @examples
##' 
##' # mean f0 one value per segment
##' m = trapply(vowlax.fund, mean, simplify=TRUE)
##' # mean F1 - F4
##' m = trapply(vowlax.fdat, apply, 2, mean, simplify=TRUE)
##' # make a logical vector of any segments that have an F1 value
##' # between their start time and end time greater than n Hz
##' pfun <- function(x, n=1000) any(x > n)
##' # greater than 1100 Hz
##' temp = trapply(vowlax.fdat[,1], pfun, 1100, simplify=TRUE)
##' # get the F2-range per segment
##' r = trapply(vowlax.fdat[,2], range)
##' # F2-range of 20th segment
##' r[[20]]
##' # DCT-smooth F2 with 10 coeffs
##' # get the first 4 DCT coefficients
##' f2.dct = trapply(vowlax.fdat[,2], dct, 3, simplify=TRUE)
##' # dct-smooth F2 with the first 5 DCT coeffs
##' f2sm = trapply(vowlax.fdat[,2], dct, 4, TRUE,  returntrack=TRUE)
##' # Make new F2 trackdata such that each segment has
##' # F2 divided by its F2 range
##' pfun <- function(x) x/(diff(abs(range(x))))
##' newf2 = trapply(vowlax.fdat[,2], pfun, returntrack=TRUE)
##' 
##' @export trapply
`trapply` <- function (trackdata, fun, ..., simplify = FALSE, returntrack = FALSE) 
{
  if(returntrack)
    simplify <- FALSE
  # if simplify is F or if returntrack is T, store as a list
  if (!simplify) 
    result <- list(NULL)
  else result <- NULL
  for (j in 1:nrow(trackdata)) {
    if (!simplify) 
      result[[j]] <- fun(trackdata[j, ]$data, ...)
    else 
      result <- rbind(result, fun(trackdata[j, ]$data, 
                                  ...))        
  }
  if (simplify) {
    if (ncol(result) == 1) 
      result <- c(result)
  }
  if (returntrack)
    result <- buildtrack(result)
  result
}
