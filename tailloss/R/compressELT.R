#' Compress the event loss table
#'
#' Function to merge losses of the same amount adding up their corresponding occurrence rates, and to round the losses to the 10^\code{digits} integer value.
#'
#' @param ELT Data frame containing two numeric columns. The column \code{Loss} contains the expected losses from each single occurrence of event. The column \code{Rate} contains the arrival rates of a single occurrence of event.
#' @param digits Integer. It specifies the rounding of the losses to the 10^\code{digits} integer value of the event loss table. \code{digits} < 0 decreases the precision of the calculation, but considerably decreases the time to perform it. If  \code{digits} = 0 it only merges the losses of the same amount adding up their corresponding rates. The default value is \code{digits} = 0.
#'  @return Data frame containg two numeric columns. The column \code{Loss} contains the expected losses from each single occurrence of event. The column \code{Rate} contains the arrival rates of a single occurrence of event.
#' @export
#' @examples
#' data(UShurricane)
#'
#' # Compress the table to thousands of dollars
#'
#' USh.k <- compressELT(ELT(UShurricane), digits = -3)
#' summary(USh.k)
#'
#' # Compress the table to millions of dollars
#'
#' USh.m <- compressELT(ELT(UShurricane), digits = -6)
#' summary(USh.m)

compressELT <- function (ELT, digits = 0) 
{
  stopifnot(inherits(ELT, "ELT"),
            length(digits) == 1, digits == floor(digits))

  ## the outer round squashes numerical noise
  
  if (!is.na(digits))
    ELT$Loss <- round(10^digits * round(ELT$Loss, digits))

  out <- aggregate(Rate ~ Loss, data = ELT, FUN = sum)
  out <- out[order(out$Loss), ]
  if (out$Loss[1] == 0)
    out <- out[-1, ] # trim off a zero
  out <- cbind(ID = 1L:nrow(out), out[, c("Rate", "Loss")])
  row.names(out) <- NULL
  attr(out, "digits") <- digits # useful to record this
  class(out) <- c("ELT", "data.frame")
  out
}
