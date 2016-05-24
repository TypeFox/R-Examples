#' @title Tajima's D
#' @description Calculate Tajima's D for a set of sequences to test 
#'   for selection.
#' 
#' @param x set of DNA sequences or a haploid \linkS4class{gtypes} 
#'   object with sequences.
#' 
#' @return A named vector with the estimate for \code{D} and 
#'   the \code{p.value} that it is different from 0.
#' 
#' @references Tajima, F. 1989. Statistical method for testing the neutral 
#'   mutation hypothesis by DNA polymorphism. Genetics 123:585-595.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom stats qbeta integrate
#' @export
#' 
tajimasD <- function(x) {
  x <- as.multidna(x)
  
  result <- do.call(rbind, lapply(getSequences(x, simplify = FALSE), function(dna) {
    dna <- as.matrix(dna)
    num.sites <- ncol(dna)
   
    pws.diff <- dist.dna(dna, model = "N", pairwise.deletion = TRUE, as.matrix = TRUE)
    pi <- mean(pws.diff[lower.tri(pws.diff)])
    S <- ncol(variableSites(dna)$site.freqs) # number segregating sites
    n <- nrow(dna) # number individuals
    
    # now for all the pieces...
    n.vec <- 1:(n-1) # common vector used
    a1 <- sum(1 / n.vec)
    a2 <- sum(1 / n.vec ^ 2)
    b1 <- (n + 1) / (3 * (n - 1))
    b2 <- 2 * (n ^ 2 + n + 3) / (9 * n * (n - 1))
    c1 <- b1 - 1 / a1
    c2 <- b2 - (n + 2) / (a1 * n) + a2 / a1 ^ 2
    e1 <- c1 / a1
    e2 <- c2 / (a1 ^ 2 + a2)
    
    # which allows you to calculate D!
    D_obs <- (pi - S / a1) / sqrt(e1 * S + e2 * S * (S - 1))
    
    if(is.nan(D_obs)) {
      warning("D cannot be computed (division by zero in final equation)")
      return(c(D = NA, p.value = NA))
    }
    
    # making life easier:
    D.to.x <- function(D, Dmin, Dmax) (D - Dmin) / (Dmax - Dmin) # conversion 1
    x.to.D <- function(x, Dmin, Dmax) x * (Dmax - Dmin) + Dmin   # conversion 2
    
    # compare with expected value from beta distribution:
    DMin <- (2 / n - 1 / a1) / e2 ^ (1 / 2) #S approaches infinity
    DMax <- if(n %% 2 == 0) {
      (n / (2 * (n - 1)) - 1 / a1) / e2 ^ (1 / 2)
    } else {
      ((n + 1)/(2 * n) - 1 / a1) / e2 ^ (1 / 2)
    }
    Alpha <- -(1 + DMin * DMax) * DMax / (DMax - DMin)
    Beta <- (1 + DMin * DMax) * DMin / (DMax - DMin)
    
    # distribution function from Tajima paper:
    beta.D <- function(tajima.D, alpha = Alpha, beta = Beta, 
                       Dmin = DMin, Dmax = DMax) { 
      # give P(D, given the other 4)
      gamma(alpha + beta) * (Dmax - tajima.D) ^ (alpha - 1) * 
        (tajima.D - Dmin) ^ (beta - 1) / 
        (gamma(alpha) * gamma(beta) * (Dmax - Dmin) ^ (alpha + beta - 1))
    }
    
    # 95% confidence interval:
    LB <- x.to.D(qbeta(.025, Beta, Alpha), DMin, DMax)
    UB <- x.to.D(qbeta(.975, Beta, Alpha), DMin, DMax)
    
    # important probabilities 
    # D negative: intregrate from Dmin to D, D positive: integrate from D to Dmax
    tryCatch({
      prob <- integrate(beta.D, lower = ifelse(D_obs < 0, DMin, D_obs),
                        upper = ifelse(D_obs < 0, D_obs, DMax), 
                        alpha = Alpha, beta = Beta, Dmin = DMin, Dmax = DMax)
      c(D = D_obs, p.value = prob$value)
    }, error = function(e) {
      warning("error in Tajima's D integration, NA returned")
      c(D = NA, p.value = NA)
    })
  }))
  
  rownames(result) <- getLocusNames(x)
  return(result)
}
