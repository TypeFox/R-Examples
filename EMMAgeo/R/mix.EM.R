#' Function to mix sample spectres.
#' 
#' This functions allows to mix grain-size distributions with specified
#' proportions and defined noise levels, for example to test the goodness of
#' the EMMA algorithm.
#' 
#' The function multiplies each end-member with the respective proportion
#' value, sums the resulting variables, adds uniform noise and normalises the
#' resulting mixed sample to 100 \%.
#' 
#' @param EM Numeric matrix containing the grain-size distribution definitions.
#' Each definition is in a separate row with variable contributions in columns.
#' @param proportion Numeric vector containing the relative proportions of each
#' distribution per sample.
#' @param noise Numeric scalar containing optional relative white noise levels.
#' @param autocorrelation Numeric scalar specifying the degree of
#' autocorrelation among classes.  Autocorrelation is realised as running mean
#' of the specified length. Only odd values are allowed.
#' @return Numeric vector comprising a sample composed of known proportions of
#' end-members.
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{create.EM}}
#' @keywords EMMA
#' @examples
#' 
#' ## define end-member loadings and phi vector
#' EMa.1 <- create.EM(p1 = c(2, 8), p2 = c(1, 0.8), s = c(0.7, 0.3), 
#'                    boundaries = c(0, 11), n = 80)
#' EMa.2 <- create.EM(p1 = c(4, 7), p2 = c(1.1, 1.4), s = c(0.5, 0.5),
#'                    boundaries = c(0, 11), n = 80)
#' EMa   <- rbind(EMa.1, EMa.2)
#' 
#' phi   <- seq(0, 11, length.out = 80)
#' 
#' ## mix end-member loadings
#' sample1 <- mix.EM(EMa, proportion = c(0.3, 0.7))
#' sample2 <- mix.EM(EMa, proportion = c(0.5, 0.5), noise = 0.1,
#'                   autocorrelation = 3)
#' 
#' ## plot end-member loadings (grey) and resulting samples (black)
#' plot(phi, EMa.1, type="l", col = "grey")
#' lines(phi, EMa.2, col = "grey")
#' lines(phi, sample1)
#' lines(phi, sample2)
#' 
#' @export mix.EM
mix.EM <- function(
  EM,
  proportion, 
  noise,
  autocorrelation
){
  
  ## check/set noise variable
  if(missing(noise) == TRUE) {noise <- 0}
  
  ## rescale individual end-member loadings according to their proportion
  EMs <- EM * as.numeric(proportion)
  
  ## sum up all individual end-member loadings to one composite
  EMs <- apply(EMs, 2, sum)
  
  ## optionally add noise
  EMs <- EMs + noise * runif(length(EMs), -max(EMs), max(EMs))
  
  ## Apply running mean to generate autocorrelation
  if(missing(autocorrelation) != TRUE) {
    if(autocorrelation %% 2 == 0) {
      stop("Value for autocorrelation is no odd integer.")
    }
    b.h <- autocorrelation / 2 - 0.5
    running.mean <- EMs
    for(i in (1 + b.h):(length(EMs) - b.h)) {
      running.mean[i] <- mean(EMs[(i - b.h):(i + b.h)])
    }
    EMs <- running.mean
  }
  
  ## set negative class contents to zero
  EMs[EMs < 0] <- 0
  
  ## rescale end-member loadings to sum up to 1
  EMs <- EMs / sum(EMs)

  return(EMs)
}