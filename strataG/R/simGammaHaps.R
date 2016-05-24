#' @title Simulate Haplotypes
#' @description Simulate a haplotypic frequency distribution based on
#'   a specified gamma distribution.
#'
#' @param pop.size size of population.
#' @param num.haps number of haplotypes to generate.
#' @param shape,scale parameters of Gamma distribution (see \code{\link{dgamma}}).
#' @param return.freq logical. Return frequency table of haplotypes? If \code{FALSE} return vector of haplotypes.
#' @param plot logical. Show plot of haplotypic frequency distribution?
#'
#' @return Frequency table of haplotypes.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' haps <- simGammaHaps(1000, 15, 1, 2.5)
#' print(haps)
#'
#' @importFrom graphics layout par curve hist box
#' @importFrom stats dgamma
#' @export

simGammaHaps <- function(pop.size, num.haps, shape, scale,
                         return.freq = TRUE, plot = TRUE) {
  if(num.haps > pop.size) stop("'num.haps' must be smaller than 'pop.size'.")

  # get frequency of each haplotype from gamma
  hap.freq <- dgamma(1:num.haps, shape = shape, scale = scale)
  hap.freq <- ceiling(hap.freq * pop.size / sum(hap.freq))
  
  # create "individuals"
  haplotypes <- rep(1:num.haps, times = hap.freq)
  haplotypes <- sample(haplotypes, pop.size)

  if(plot) {
    gamma.func <- function(x) dgamma(x, shape = shape, scale = scale)
    layout(matrix(c(1, 2), nrow = 2))
    op <- par(mar = c(4, 4, 3, 1) + 0.1)
    curve(gamma.func, 0, num.haps,
          xlab = "Haplotype", ylab = "Density",
          main = paste("Gamma(", shape, ", ", scale, ")", sep = "")
    )
    par(mar = c(5, 4, 1, 1) + 0.1)
    hist(haplotypes, xlab = "Haplotype", ylab = "Density", main = "")
    box()
    layout(matrix(1, 1))
    par(op)
  }

  if(return.freq) table(haplotypes) else haplotypes
}
