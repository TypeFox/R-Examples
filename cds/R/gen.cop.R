#' Generate a Copula
#' 
#' Generate correlated data multivariate categorical data via a copula. 
#' 
#' @param n Integer; the number of samples to draw.
#' @param tauvek A vector of association parameters for each of the Clayton copulae 
#' (see \code{\link{copClayton}}), of the same length as \code{nr.cols}.
#' @param nr.cols A vector giving the number of columns to draw from each of the copulae.
#' @param true.mu A vector giving the mean for each of the columns in the data.
#' @param err.coeff The standard errors for underlying normal distribution.
#' @param random Logical indicating whether or not the samples should be presented in 
#' random order.
#' @param reverse Logical indicating whether some of the simulated variables should be reversed to 
#' have negative association or not.
#' @param reverse.thresh The proportion of columns to reverse.
#' @keywords multivariate
#' @export gen.cop
gen.cop <-
function(n, tauvek = c(0.2, 0.35), nr.cols = c(10, 10), 
                    true.mu = runif(sum(nr.cols)), err.coeff = 0.1, random = FALSE, 
                    reverse = TRUE, reverse.thresh = 0.75)
{
  env.id <- environment()
  m <- sum(nr.cols)
  if(length(tauvek) != length(nr.cols)) stop("Length of tauvek and nr.cols does not match")
  nr <- length(tauvek)
  csum <- c(0, cumsum(nr.cols))
  samp <- matrix(NA, ncol = m, nrow = n)
  true.tau <- matrix(0, nrow = m, ncol = m)
  for (i in seq_len(nr)) {
    tmp <- substitute(copula::onacopula("Clayton", C(copula::copClayton@iTau(tauvek[i]), 1:nr.cols[i])))
    Cop <- eval(tmp)
    samp[,(csum[i]+1):csum[i+1]] <- copula::rnacopula(n, Cop)
    true.tau[(csum[i]+1):csum[i+1],(csum[i]+1):csum[i+1]] <- tauvek[i]
  }
  diag(true.tau) <- 1
  if (random) {
    ord <- sample(1:ncol(samp), ncol(samp))
    samp <- samp[,ord]
    true.tau <- true.tau[ord,ord]
  }
  if (reverse) {
    ind <- runif(sum(nr.cols)) > reverse.thresh
    samp[, ind] <- 1 - samp[,ind]
    true.tau[ind, ] <- -1*true.tau[ind,]
    true.tau[-ind, ind] <- -1*true.tau[-ind,ind]
    diag(true.tau) <- 1
  } 
  samp <- apply(rbind(true.mu, samp), 2, function(x) trQnorm(x[-1], 
                mean = x[1], sd = err.coeff))
  list(sample = samp, tau = true.tau)
}
