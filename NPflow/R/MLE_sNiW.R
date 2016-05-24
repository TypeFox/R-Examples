#'MLE for sNiW distributed observations
#'
#'Maximum likelihood estimation of Normal inverse Wishart distributed observations
#'
#'@rdname MLE_sNiW
#'
#'@param xi_list a list of length \code{N} whose elements are observed vectors of length \code{d}
#'of the mean parameters xi.
#'
#'@param psi_list a list of length \code{N} whose elements are observed vectors of length \code{d}
#'of the skew parameters psi.
#'
#'@param S_list a list of length \code{N} whose elements are observed variance-covariance matrices
#'of dimension \code{d x d}.
#'
#'@param doPlot a logical flag indicating wether the algorithm progression should be plotted.
#'Default is \code{TRUE}.
#'
#'@importFrom stats uniroot
#'
#'@author Boris Hejblum, Chariff Alkhassim
#'
#'@export
#'
#'@examples
#'hyperG0 <- list()
#'hyperG0$b_xi <- c(0.3, -1.5)
#'hyperG0$b_psi <- c(0, 0)
#'hyperG0$kappa <- 0.001
#'hyperG0$D_xi <- 100
#'hyperG0$D_psi <- 100
#'hyperG0$nu <- 35
#'hyperG0$lambda <- diag(c(0.25,0.35))
#'
#'xi_list <- list()
#'psi_list <- list()
#'S_list <- list()
#'for(k in 1:1000){
#'  NNiW <- rNNiW(hyperG0, diagVar=FALSE)
#'  xi_list[[k]] <- NNiW[["xi"]]
#'  psi_list[[k]] <- NNiW[["psi"]]
#'  S_list[[k]] <- NNiW[["S"]]
#'}
#'
#'mle <- MLE_sNiW(xi_list, psi_list, S_list)
#'mle

MLE_sNiW <- function(xi_list, psi_list, S_list, doPlot=TRUE){


  N <- length(xi_list)
  d <- length(xi_list[[1]])

  if(length(psi_list) != N | length(S_list) != N){
    stop("Number of observations/MCMC iterations is not matching")
  }

  Sinv_list <- lapply(S_list, solve)
  Sinv_sum <- Reduce('+', Sinv_list)

  U_xi_U_psi<-Reduce('+',lapply(1:N,
                                function(i,m) {rbind(xi_list[[i]],psi_list[[i]])%*%m[[i]]},
                                m=Sinv_list))%*% solve(Sinv_sum)
  U_xi<-U_xi_U_psi[1,]
  U_psi<-U_xi_U_psi[2,]

  xim <- lapply(xi_list, function(x){x - U_xi})
  psim <- lapply(psi_list, function(x){x - U_psi})

  U_B <- N*d*solve(matrix(rowSums(mapply(x = xim, p = psim, Si = Sinv_list, FUN=function(x,p,Si){
    v <- rbind(x, p)
    tcrossprod(v%*%Si,v)
  }, SIMPLIFY=TRUE)),
  nrow=2, byrow=FALSE))

  U_df<- try(stats::uniroot(function(nu0){(N/2*digamma_mv(x=nu0/2, p=d)
                                           + 1/2*sum(sapply(S_list, function(S){log(det(S))}))
                                           - N*d/2*log(N*nu0/2)
                                           + N/2*log(det(Sinv_sum))
  )}, lower = d+1, upper=1E9)$root, TRUE)
  if(inherits(U_df, "try-error")){U_df <- d+1}

  U_Sigma <- N*U_df*solve(Sinv_sum)


  return(list("U_xi" = U_xi,
              "U_psi" = U_psi,
              "U_B" = U_B,
              "U_df" = U_df,
              "U_Sigma" = U_Sigma))

}