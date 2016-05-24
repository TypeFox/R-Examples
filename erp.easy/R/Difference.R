#' Calculate a difference waveform from two data frames
#'
#' \code{dif.wave} calculates a difference waveform from two data frames in the format
#'   returned from \code{load.data}
#'
#' @param x A data frame in the format returned from \code{\link{load.data}}. This value serves
#'   as the minuend (i.e., value to be subtracted from).
#' @param z A data frame in the format returned from \code{load.data}. This value serves
#'   as the subtrahend (i.e., value subtracted).
#' @param name Specify the Stimulus column name of the difference data frame.  Must provide the new
#'   name in quotations.
#' @param keep Include one of the data frames used in subtraction in the returned difference data frame.
#'   By default \code{keep} = NULL and only the resulting difference data frame will be returned.
#'
#' @details The data frames must either: \enumerate{
#'    \item have been imported using \code{load.data}
#'    - OR -
#'    \item for subset data, the "Stimulus" column must be "refactored." See the example
#'    below for one way to "reset" the factors in the "Stimulus" column.  Note this procedure
#'    is only necessary for subset data, and not data frames that were independently imported.
#'}
#'
#' @return One of two possible data frames, depending on the value of \code{keep}:
#' \enumerate{
#'   \item \code{keep} = "y": the difference (i.e., x - z) data frame
#'     and the minuend data frame (i.e., x)
#'   \item \code{keep} = "n": just the difference data frame
#'}
#'
#' @examples
#'  # Calculate a difference wave
#'  Negative = ERPdata[1:6765, ]
#'  Neutral = ERPdata[6766:13530, ]
#'  refactor.neg <- factor(Negative$Stimulus)
#'  refactor.neut <- factor(Neutral$Stimulus)
#'  Negative$Stimulus <- refactor.neg
#'  Neutral$Stimulus <- refactor.neut
#'  difference <- dif.wave(Negative, Neutral, name = "Neg - Neut", keep = Neutral)
#'  head(difference) # view the new Stimulus column name
#'  grandaverage(difference, "V78") # plot the grand average difference wave
#'
#' @author Travis Moore

dif.wave <- function(x, z, name = NULL, keep = NULL) {
  fcolsx <- x[ , 1:3]
  fcolsz <- z[ , 1:3]
  xcols <- x[ , 4:ncol(x)]
  zcols <- z[ , 4:ncol(z)]
  diff.cols <- xcols - zcols
  final <- cbind.data.frame(fcolsz, diff.cols)

  # resets the Stimulus column factors
    if (!is.null(name)) {
      final$Stimulus <- name
    } else {
      new.name <- paste(x[1,2], " - ", z[1,2], sep = "")
      final$Stimulus <- new.name
    }
  refactor <- factor(final$Stimulus)
  final$Stimulus <- refactor
  if (!is.null(keep)) {
    combine <- rbind(final, keep)
    return(combine)
  } else {
    return(final)
  }
}
