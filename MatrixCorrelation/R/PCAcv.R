#' @title Principal Component Analysis cross-validation error
#'
#' @description PRESS values for PCA as implemented by Eigenvector and described by Bro et al. (2008).
#'
#' @param X \code{matrix} object to perform PCA on.
#' @param ncomp \code{integer} number of components.
#'
#' @details For each number of components predicted residual sum of squares are calculated
#' based on leave-one-out cross-validation. The implementation ensures no over-fitting or
#' information bleeding.
#'
#' @return A vector of PRESS-values.
#'
#' @author Kristian Hovde Liland
#'
#' @references R. Bro, K. Kjeldahl, A.K. Smilde, H.A.L. Kiers, Cross-validation of component models: A critical look at current methods. Anal Bioanal Chem (2008) 390: 1241-1251.
#'
#' @seealso \code{\link{plot.SMI}} (print.SMI/summary.SMI), \code{\link{RV}} (RV2/RVadj), \code{\link{r1}} (r2/r3/r4/GCD), \code{\link{allCorrelations}} (matrix correlation comparison).
#'
#' @examples
#' X1  <- scale( matrix( rnorm(100*300), 100,300), scale = FALSE)
#' PCAcv(X1,10)
#'
#' @export
PCAcv <- function(X, ncomp){

  X <- as.matrix(X)
  class(X) <- "matrix"
  N <- dim(X)
  
  # Center X
  X <- X - rep(colMeans(X), each = N[1])
  
  # Set ncomp
  if(missing(ncomp)){
    ncomp <- min(N[1]-1,N[2])
  } else {
    ncomp <- min(ncomp,min(N[1]-1,N[2]))
  }
  
  # Prepare storage
  Xhat <- array(0, dim = c(N[1],N[2],ncomp))
  
  # Cross-validation (leave-one-out)
  pb <- progress_bar$new(total = N[1], format = "  [:bar] :percent (:eta)")
  for(i in 1:N[1]){
    Xi  <- X[-i, , drop = FALSE]
    Pi  <- svds(Xi, k = ncomp,  nv = ncomp, nu = 0)$v
    Xii <- matrix(rep(X[i,], N[2]), N[2], N[2], byrow = TRUE)
    diag(Xii) <- 0
    
    # Magic to avoid information bleed
    PiP  <- apply(Pi^2, 1, cumsum)
    PiP1 <- t(PiP/(1-PiP)+1)
    PihP <- t(Pi*(Xii%*%Pi))
    for(j in 1:N[2]){
      PP <- PihP[,j, drop = FALSE] %*% PiP1[j,, drop = FALSE]
      PP[lower.tri(PP)] <- 0
      Xhat[i,j, ] <- colSums(PP)
    }
    pb$tick()
  }
  
  error <- numeric(ncomp)
  for(i in 1:ncomp){
    error[i] <- sum((X-Xhat[,,i])^2)
  }
  error
}
