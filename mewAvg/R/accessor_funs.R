#' @title Extract MEW average value
#'
#' @description Return the current value of the moving expanding
#' window (MEW) average if it is up-to-date; otherwise, raise an error
#'
#' @param av The current state of the MEW average
#'
#' @return (vector double length n_xx) the current value of the MEW
#' average if it is up-to-date
#'
#' @examples
#' ## see the examples for the function \code{mewMean}
#'
#' @export
mewGetMean <- function (av) {

  if (class(av)[1] != "mewTyp") {

    stop("mewGetMean: argument must have class mewTyp")
  }

  if (av@know_mean == 0) {

    stop("mewGetMean: the MEW average is NOT up-to-date")
  } else if (av@know_mean == 1){

    return(av@x_mean)
  } else {

    stop("mewGetMean: know_mean is set inappropriately")
  }
}
