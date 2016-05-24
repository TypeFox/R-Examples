#' Grouped Simulation with Response Styles
#' 
#' Simulate response data for a group of response styles.
#' 
#' @param nr.indv A vector giving the number of respondents in each group.
#' @param m The number of objects.
#' @param scales The rating scale used, 1:q.
#' @param err.coeff The standard error used in the underlying normal noise.
#' @param alphamat The matrix of spline parameters defining the response styles, with each
#' row containing a response style. No intercepts should be included.
#' @param true.mu Optional; a matrix or vector giving the true underlying preferences for the objects.
#' @param random Logical indicating whether to apply the response styles in random order
#' @param same.mu Logical indicating whether a universal value for mu should be assumed.
#' @param use.copula Logical indicating whether to use a correlated dependence structure 
#' through a copula.
#' @param reverse.thresh A numeric value giving the proportion of observations for which the 
#' dependece structure should be reversed. Only applicable when \code{copula} is \code{TRUE}.
#' @return An object of class \code{cdsdata}, inheriting from class \code{icdsdata}, which is a 
#' list with the following slots:
#' \describe{
#'  \item{prers}{The pre-response style simulated data}
#'  \item{postrs}{The data after adding the response styles}
#'  \item{postbl}{The same as \code{postrs} in this case}
#'  \item{Fr.cent.rs}{The centred Fr matrix for \code{postrs}}
#'  \item{Fr.rs}{The Fr matrix for \code{postrs}}
#'  \item{Fr.cent.bl}{The same as \code{Fr.cent.rs}, for compatibility with \code{icds}}
#'  \item{Fr.bl}{The same as \code{Fr.rs}, for compatibility with \code{icds}}
#'  \item{mu}{Matrix of the true underlying preference structure for the obects}
#'  \item{block}{Numeric vector identifying the different blocks for incompleteness, in this case a vector
#'  of ones}
#'  \item{grp.rs}{The response style grouping vector}
#'  \item{alphamat}{Matrix of spline parameters for the response styles}
#'  \item{scales}{The rating scale 1:q used}
#'  \item{m}{Number of objects}
#'  \item{munique}{The number of objects seen within each block - equal to zero in this case}
#'  \item{m0}{The number of objects seen by all subjects - equal to \code{m} in this case}
#'  \item{true.tau}{Actual tau used in the simulation with copulae}
#'  \item{call}{The function call}
#' }
#' 
#' @seealso \code{\link{createcdsdata}}
#' @keywords multivariate
#' @export cds.sim
cds.sim <- function(nr.indv = c(100, 100, 100), m = 25, scales = 1:7, err.coeff = 0.1, 
         alphamat = rbind(c(4, 4, 1), c(1, 4, 4), c(1, 2, 1)), true.mu = NULL, 
         random = TRUE, same.mu = TRUE, use.copula = FALSE,
         reverse.thresh = 1)
{
  cll <- match.call()
  if(use.copula) {
    cat("Using copula: number of objects m set to 20\n")
    m <- 20
  }
  K <- length(nr.indv)
  n <- sum(nr.indv)
  out.lst <- vector("list", length = K)
  q <- length(scales)
  mumat <- matrix(NA, nrow = K, ncol = m)
  if(same.mu) true.mu <- runif(m)
  for(k in 1:K) {
    out <- datsim(resp.style = create.rs(alpha = alphamat[k,], nr.scale = q), 
                  nr.indv = nr.indv[k], m = m, scales = scales, 
                  err.coeff = err.coeff, true.mu = true.mu, 
                  use.copula = use.copula, reverse.thresh = reverse.thresh)
    
    out.lst[[k]] <- out
    if(k == 1) {
        Tmat <- out$Tmat
        prers <- out$prers
        postrs <- out$postrs
      } else {
        Tmat <- rbind(Tmat, out$Tmat)
        prers <- rbind(prers, out$prers)
        postrs <- rbind(postrs, out$postrs)
      }
    mumat[k,] <- out$mu
  }
  colnames(prers) <- colnames(postrs) <- paste0("Item", 1:m)
  
  ## Randomly reorder
  grp.rs <- rep(1:K, nr.indv)
  if(random) {
        ord <- sample(1:n,  n)
        grp.rs <- grp.rs[ord]
        Tmat <- Tmat[ord,]
        postrs <- postrs[ord,]
        prers <- prers[ord,]
      }
  
  ## Calculate Fr and Fr.cent
  Fr <- rbind(Tmat, m + q - 2 - Tmat)
  Fr.cent <- Fr - 0.5*(m + q - 2) * matrix(1, nrow = 2*n, ncol = m + q - 1)
  
  colnames(alphamat) <- paste0("a", 1:3)
  out <- list(prers = prers, postrs = postrs, postbl = postrs, Fr.cent.rs = Fr.cent, Fr.rs = Fr, 
              Fr.cent.bl = Fr.cent, Fr.bl = Fr, block = rep(1, n), 
              mu = mumat, grp.rs = grp.rs, alphamat = alphamat,  
              scales = scales,  m = m, m0 = m, munique = 0, true.tau = out$tau, call = cll)
  class(out) <- c("cdsdata",  "icdsdata", "list")
  out
}
