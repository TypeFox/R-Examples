#' Simple and stratified systematic sampling
#' @description Select sampling units using simple or stratified systematic samplin. In the context of two-stage cluster sampling, select Secondary Sampling Units (SSU) in one or more Primary Sampling Units (PSU), using systematic sampling.
#' @param psu.ssu \code{\link{data.frame}} with all PSU. First column contains PSU unique identifiers. Second column contains \code{\link{numeric}} PSU sizes. Only used for the second stage of a two-stage cluster design (see details).
#' @param su \code{\link{numeric}} indicating the number of sampling units to be selected. If \code{su} has more than one element, stratified sampling is applied and \code{psu.ssu} is ignored (see details).
#' @param N \code{\link{numeric}} indicating the number of sampling units in the population. It is intended for simple or stratified sampling designs and when used, \code{psu.ssu} is ignored (see details).
#' @param write logical. If \code{TRUE}, a *.csv file containing the PSU and their SSU is writed in the current working directory.
#' @param ... further arguments passed to \code{\link{write.table}} function.
#' @return A \code{matrix}. For the second stage in a two-stage cluster sampling, the names of columns are the identifiers of selected psu, coerced by \code{\link{as.character}} to avoid scientific notation in case the identifiers be of \code{\link{class}} \code{\link{numeric}}. The rows correspond to the selected SSU within each psu. For simple systematic sampling, the rows correspond to the selected sampling units. For stratified sampling, each column represent an strata and the rows correspond to the selected sampling units in each strata.
#' @details When \code{N} is defined, \code{psu.ssu} is ignored (not need to be defined). If \code{N} has one element, \code{su} must too and the result is a simple systematic selection.  If \code{N} has more than one element, \code{su} must have the same number of elements and each oredered pair represent an strata. Thus, when N has more than one element, the result is a stratified sampling with systematic selection within each strata (see examples).
#' @references Levy P and Lemeshow S (2008). Sampling of populations: methods and applications, Fourth edition. John Wiley and Sons, Inc.
#' 
#' \url{http://oswaldosantos.github.io/capm}
#' @seealso \code{\link{SamplePPS}}.
#' @export
#' @examples 
#' # Load data with PSU identifiers and sizes.
#' data(psu.ssu)
#' 
#' # Take a sample of 10 PSU, with probability 
#' # proportional to size and with replacement.
#' selected.psu <- SamplePPS(psu.ssu, 10, write = FALSE)
#' 
#' # Take a systematic sampling of 5 SSU within each 
#' # PSU of selected.psu.
#' SampleSystematic(selected.psu, 5, write = FALSE)
#' 
#' # Simple systematic sampling
#' SampleSystematic(su = 5, N = 100)
#' 
#' # Stratified systematic sampling
#' SampleSystematic(su = c('Urban' = 50, 'Rural' = 10),
#'                  N = c('Urban' = 4000, 'Rural' = 150))
SampleSystematic <- function(psu.ssu = NULL, su = NULL, N = NULL, write = FALSE, ...) {
  if (!is.null(N)) {
    sus <- matrix(rep(NA, max(su) * length(su)), ncol = length(N))
    for (i in 1:length(N)) {
      if (N[i] <= su[i]) {
        stop('Sampling units in the sample must be lesser than in population.')
      }
      if (length(N) != length(su)) {
        stop('N and su must have the same number of elements.')
      }
      int <- floor(N[i] / su[i])
      k <- sample(int, 1)
      Su <- rep(NA, su[i])
      for(j in 1:su[i]) {
        Su[j] <- k
        k <- k + int
      }
      sus[ , i][1:su[i]] <- floor(Su)
    }
    if (!is.null(names(N))) {
      colnames(sus) <- as.character(names(N)) 
    }
    if (write) {
      write.table(sus, file = 'sampling_units.csv', sep = ',', dec = '.', 
                  qmethod = 'double', col.names = NA, ...)
      cat('\n', 'The \"sampling_units.csv\" file contains the selected. \nsampling units. It is in the directory:', '\n\n', getwd(), '\n', '\n')
    }
    return(sus)
  } else {
    inv <- c(which(!is.finite(psu.ssu[, 2])), which(psu.ssu[, 2] <= 0))
    if (length(inv) > 0) {
      stop('The size of the following sampling unit(s) is(are) invalid:', '\n', paste('   ', inv))
    }
    inv1 <- which(psu.ssu[, 2] < su)
    if (length(inv1) > 1) {
      stop('The number of secondary sampling units to be selected (', su, ') is greater than the size of the following primary sampling units:', '\n', paste('   ', inv1))
    }
    sus <- matrix(rep(NA, nrow(psu.ssu) * su), ncol = nrow(psu.ssu))
    int <- floor(psu.ssu[, 2] / su)
    for(i in 1:length(int)) {
      k <- sample(int[i], 1)
      Su <- rep(NA, su)
      for(j in 1:su) {
        Su[j] <- k
        k <- k + int[i]
      }
      sus[ , i] <- floor(Su)
    }
    colnames(sus) <- as.character(psu.ssu[ , 1])
    if (write == TRUE) {
      colnames(sus) <- as.character(psu.ssu[ , 1])
      write.table(sus, file = 'sus.csv', sep = ',', dec = '.', 
                  qmethod = 'double', col.names = NA, ...)
      cat('\n', 'The \"sus.csv\" file contains the selected', '\n', 'sampling units. It is in the directory:', '\n\n', getwd(), '\n', '\n')
    }
    return(sus)
  }
}
