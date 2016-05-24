#' Probability proportional to size sampling with replacement
#' @description Select Primary Sampling Units (PSU) with probability proportional to size with replacement.
#' @param psu.ssu \code{\link{data.frame}} with all PSU. First column contains PSU unique identifiers. Second column contains \code{\link{numeric}} PSU sizes.
#' @param psu the number of PSU to be selected.
#' @param write logical. If \code{TRUE}, a *.csv file containing the PSU and their Secondary Sampling Units (SSU) is writed in the current working directory.
#' @param ... further arguments passed to \code{\link{write.table}} function.
#' @return \code{\link{data.frame}}. First column contains the selected psu identifiers, coerced by \code{\link{as.character}}, to avoid scientific notation in case the identifiers be large numbers of \code{\link{class}} \code{\link{numeric}}. Second column contain PSU sizes, a variable needed for second stage sampling with \code{\link{SampleSystematic}}.
#' @references Levy P and Lemeshow S (2008). Sampling of populations: methods and applications, Fourth edition. John Wiley and Sons, Inc.
#' 
#' \url{http://oswaldosantos.github.io/capm}
#' @seealso \code{\link{SampleSystematic}}.
#' @export
#' @examples 
#' # Load data with PSU identifiers and sizes.
#' data(psu.ssu)
#' 
#' # Take a sample of 10 PSU with probability proportional to size with replacement.
#' (selected.psu <- SamplePPS(psu.ssu, 10, write = FALSE))
#' 
SamplePPS <- function (psu.ssu = NULL, psu = NULL, write = FALSE, ...) {
  inv <- c(which(!is.finite(psu.ssu[, 2])), which(psu.ssu[, 2] <= 0))
  if (length(inv) > 0) {
    stop('The size of the following sampling unit(s) is(are) invalid:', '\n', paste('   ', inv))
  }
  if (psu < 1) {
    return(NULL)
  }
  if (psu > nrow(psu.ssu)) {
    stop('The number of sampling units to be selected (', psu, ') is greater than the total number of sampling units in the population (', nrow(psu.ssu), ').')
  }
  inv2 <- which(psu.ssu[, 1] == psu.ssu[, 1][duplicated(psu.ssu[, 1])])
  if (length(inv2) > 1) {
    stop('The following psu are repeated:', '\n', paste('   ', psu.ssu[, 1][duplicated(psu.ssu[, 1])]), '\n', 'It appears in positions:', '\n', paste('   ', inv2))
  }
  M <- nrow(psu.ssu) 
  cum <- cumsum(psu.ssu[ , 2]) 
  N <- cum[M] 
  Psu <- data.frame('selected psu' = rep(NA, psu), size = rep(NA, psu)) 
  for (i in 1:psu) { 
    a <- runif(1, 0, N) 
    j <- 1
    while (cum[j] < a) { 
      j <- j + 1
    }
    Psu[i, ] <- psu.ssu[j, ] 
  }
  Psu[, 1] <- as.character(Psu[,1])
  if (write == T) {
    write.table(Psu, file = 'selected_psu.csv', sep = ',', dec = '.', 
                qmethod = 'double', row.names = FALSE, ...)
    cat('\n', 'The \"selected_psu.csv\" file contains the selected', '\n', 'PSU and their sizes (SSU). It is in the directory:', '\n\n', getwd(), '\n', '\n')
  }
  return(Psu)
}
