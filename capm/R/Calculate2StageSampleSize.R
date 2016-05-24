#' Two-stage cluster sampling size and composition
#' @description Calculates sample size and composition for a two-stage cluster sampling design to estimate a total.
#' @param psu.ssu \code{\link{data.frame}} with all primary sampling units (PSU). First column contains PSU unique identifiers. Second column contains \code{\link{numeric}} PSU sizes.
#' @param psu.x \code{\link{data.frame}}. Each row corresponds to a secondary sampling unit (SSU) surveyed in a pilot study. First column contains the PSU identifiers to which the ssu belongs to. Second column contains the totals observed in the ssu and must be \code{\link{numeric}}.
#' @param conf.level the confidence level required. It must be \code{\link{numeric}} between 0 and 1 inclusive.
#' @param error the maximum relative difference between the estimate and the unknown population value. It must be \code{\link{numeric}} between 0 and 1 inclusive.
#' @param cost the ratio of the cost of sampling a PSU to the cost of sampling a SSU.
#' @param minimum.ssu integer to define the minimum number of SSU to be selected per PSU. If the calculated number of SSU to be selected is lesser than \code{minimum.ssu}, it is redefined as \code{minimum.ssu}. To avoid any lower threshold, define \code{minimum.ssu} as equal to 0.
#' @return Matrix with the sample size and composition and with variability estimates.
#' @details It is assumed that psu from the pilot are selected with probability proportional to size (PPS) and with replacement. ssu are assumed to be selected via simple (systematic) random sampling.
#' 
#' PSU must have the same identifiers in \code{psu.ssu} and in \code{psu.x}.
#' @references Levy P and Lemeshow S (2008). Sampling of populations: methods and applications, Fourth edition. John Wiley and Sons, Inc.
#' 
#' \url{http://oswaldosantos.github.io/capm}
#' @export
#' @examples 
#' # Load data with psu identifiers and sizes.
#' data(psu.ssu)
#' 
#' # Load data from a pilot sample.
#' data(pilot)
#' 
#' # Calculate sample size and composition.
#' (sample.sc <- Calculate2StageSampleSize(psu.ssu, pilot, conf.level = 0.95, error = 0.1, cost = 4))

Calculate2StageSampleSize <- function(psu.ssu = NULL, psu.x = NULL, conf.level = .95, error = 0.1, cost = 4, minimum.ssu = 15) {
  if (length(intersect(psu.ssu[, 1], psu.x[, 1])) == 0) {
    stop('PSU identifiers must be equal in psu.ssu and in psu.x')
  }
  if (conf.level > 1 | conf.level < 0) {
    stop('conf.level must be a number between 0 and 1 inclusive.')
  }
  if (error > 1 | error < 0) {
    stop('error must be a number between 0 and 1 inclusive.')
  }
  psu.ssu.x <- merge(psu.ssu, psu.x, by = 1)
  M <- nrow(psu.ssu)
  N <- sum(psu.ssu[ , 2])
  Ni <- psu.ssu[ , 2]
  Nip <- tapply(psu.ssu.x[, 2], psu.ssu.x[ , 1], unique)
  Nb <- mean(Ni) 
  nip <- tapply(psu.ssu.x[, 2], psu.ssu.x[ , 1], length)
  nbp <- mean(nip)
  mp <- length(unique(psu.ssu.x[ , 1])) 
  np <- nrow(psu.ssu.x)
  xi <- tapply(psu.ssu.x[ , 3], psu.ssu.x[ , 1], sum)
  Xi <- xi * Nip / nip
  vec <- sum((Xi - mean(Xi)) ^ 2) / mp
  dq <- as.numeric(unlist(tapply(
    psu.ssu.x[ , 3], psu.ssu.x[ , 1],
    function(x) (x - mean(x))^2))) 
  vdc <- sum((Nip / (Nip - 1)) * 
               (tapply(dq, psu.ssu.x[ , 1], sum))) / sum(nip)
  d <- (((M / (M - 1)) * vec) - (Nb * vdc)) / 
    (((M / (M - 1)) * vec) + (Nb * (Nb - 1) * vdc)) 
  d <- if (d < 0 | d == 0) {d = 1e-03} else {d = d}
  nb <- ceiling(sqrt(cost * ((1 - d) / d)))
  if (nb < minimum.ssu) {nb = minimum.ssu}
  X <- sum(N / sum(nip) * tapply(psu.ssu.x[ , 3], 
                                 psu.ssu.x[ , 1], sum)) 
  z <- abs(round(qnorm((1 - conf.level) / 2, 0, 1), 2))
  m <- ceiling(((z ^ 2) * sum((((N * xi) / nbp) - X) ^ 2)) / 
                 ((error ^ 2) * (X ^ 2) * (mp - 1)))
  if(m > M) {m  <-  M}
  
  sam <- matrix(c(m * nb, m, nb, vec, vdc, d), 
                ncol = 1)
  rownames(sam) <- c('Sample size',
                     'Number of PSU to be sampled',
                     'Number of SSU to be sampled in each psu',
                     'Intercluster variance',
                     'Intracluster variance',
                     'Intraclass correlation coefficient')
  colnames(sam) <- 'Value'
  return(sam)
}