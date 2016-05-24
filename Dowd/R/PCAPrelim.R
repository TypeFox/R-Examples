#' 
#' Estimates VaR plot using principal components analysis
#' 
#' @param Ra Matrix return data set where each row is interpreted as a set of 
#' daily observations, and each column as the returns to each position in a 
#' portfolio position
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#'
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Computes PCA Prelim
#'    # This code was based on Dowd's code and similar to Dowd's code,
#'    # it is inconsistent for non-scalar data (Ra).
#'    library(MASS)
#'    Ra <- .15
#'    PCAPrelim(Ra)
#'
#' @import MASS
#'
#' @export
PCAPrelim <- function(Ra){
  rho <- as.matrix(Ra)
  corr.matrix <- rbind(cbind(rho%^%0, rho%^%1, rho%^%2, rho%^%3, rho%^%4, 
                             rho%^%5, rho%^%6, rho%^%7, rho%^%8, rho%^%9),
                       cbind(rho%^%1, rho%^%0, rho%^%1, rho%^%2, rho%^%3,
                             rho%^%4, rho%^%5, rho%^%6, rho%^%7, rho%^%8),
                       cbind(rho%^%2, rho%^%1, rho%^%0, rho%^%1, rho%^%2,
                             rho%^%3, rho%^%4, rho%^%5, rho%^%6, rho%^%7),
                       cbind(rho%^%3, rho%^%2, rho%^%1, rho%^%0, rho%^%1,
                             rho%^%2, rho%^%3, rho%^%4, rho%^%5, rho%^%6),
                       cbind(rho%^%4, rho%^%3, rho%^%2, rho%^%1, rho%^%0,
                             rho%^%1, rho%^%2, rho%^%3, rho%^%4, rho%^%5),
                       cbind(rho%^%5, rho%^%4, rho%^%3, rho%^%2, rho%^%1,
                             rho%^%0, rho%^%1, rho%^%2, rho%^%3, rho%^%4),
                       cbind(rho%^%6, rho%^%5, rho%^%4, rho%^%3, rho%^%2,
                             rho%^%1, rho%^%0, rho%^%1, rho%^%2, rho%^%3),
                       cbind(rho%^%7, rho%^%6, rho%^%5, rho%^%4, rho%^%3,
                             rho%^%2, rho%^%1, rho%^%0, rho%^%1, rho%^%2),
                       cbind(rho%^%8, rho%^%7, rho%^%6, rho%^%5, rho%^%4,
                             rho%^%3, rho%^%2, rho%^%1, rho%^%0, rho%^%1),
                       cbind(rho%^%9, rho%^%8, rho%^%7, rho%^%6, rho%^%5,
                             rho%^%4, rho%^%3, rho%^%2, rho%^%1, rho%^%0))
  sigma <- corr.matrix
  mu <- double(dim(sigma)[1])
  # Random number generation
  returns <- mvrnorm(1000, mu, sigma)
  # Dowd code uses princomp in matlab. Similar function "princomp" is available 
  # in "stats" package. However, the return values from princomp are not used
  # explicitly. So, following alternative was used.
  variances <- eigen(cov(returns))$values # eigenvalues of covariance matrix.
  
  # Scree Plot
  n <- 1000
  par(c(2,1))
  percent.explained <- 100 * variances / sum(variances)
  barplot(percent.explained, xlab = "%")
  
  cum.variance <- double(length(variances))
  cum.variance[1] <- percent.explained[1]
  for (i in 2:length(variances)) {
    cum.variance[i] <- percent.explained[i] + cum.variance[i-1]
  }
  t <- 0:10
  plot(t, c(0, cum.variance), xlab = "Principal component", ylab = "%", type="l")
  title("Explanatory Power of the Principal Components")
  
}

# ------------------------------ Helper function ------------------------------
# Matrix exponentiation
"%^%" <- function(S, power) {
  # Uses eigenvalue decomposition A = PDP^-1 for matrix exponentiation.
  # Also see expm package for larger matrices for efficiency.
  y <- with(eigen(S), vectors %*% (values^power * solve(vectors)))
  return(y)
}