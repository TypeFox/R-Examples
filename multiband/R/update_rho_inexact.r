#' Update Lipschitz constant for phase majorization
#' 
#' \code{update_Lipschitz} computes a Lipschitz constant for the gradient of the objective function with the amplitude parameters fixed.
#' 
#' @param tms list of matrices whose rows are the triple (t,m,sigma) for each band.
#' @param beta intercept estimates
#' @param a amplitude estimates
#' @param inflate factor by which to multiply the Lipschitz constants
#' @export
#' @examples
#' test_data <- synthetic_multiband()
#' tms <- test_data$tms
#' beta <- test_data$beta
#' a <- test_data$a
#' update_Lipschitz(tms,beta,a)
update_Lipschitz <- function(tms,beta,a,inflate=1) {
  B <- length(tms)
  L <- double(B)
  for (b in 1:B) {
    band <- tms[[b]]
    n <- nrow(band)
    w <- 1/(band[,3]**2)
    kappa <- sum(w)
    L[b] <- abs(a[b])*abs((a[b]*kappa + sqrt(n)*norm(as.matrix((band[,2] - beta[b])/(band[,3]**2)),'f')))
  }
  return(L)
}


#' Update working response in rho update
#' 
#' \code{update_zeta} computes the working response for the MM algorithm implemented in \code{update_rho}.
#' 
#' @param tms list of matrices whose rows are the triple (t,m,sigma) for each band.
#' @param beta intercept estimates
#' @param a amplitude estimates
#' @param rho vector of the current estimates of the phase
#' @param L vector of Lipschitz constants
#' @param omega frequency
#' @export
update_zeta <- function(tms,beta,a,rho,L,omega) {
  B <- length(tms)
  zeta <- double(B)
  for (b in 1:B) {
    band <- tms[[b]]
    n <- nrow(band)
    theta <- omega*band[,1] + rho[b]
    s <- sin(theta)
    c <- cos(theta)
    m <- band[,2]
    r <- beta[b] + a[b]*s - m
    fp <- a[b]*(t(r)%*%(c/(band[,3]**2)))
    zeta[b] <- L[b]*rho[b] - fp
  }
  return(zeta)
}

#' Majorization for phase
#' 
#' \code{mm_phase_obj} computes the convex quadratic majorization used in the phase estimation step.
#' 
#' @param rho vector of new phase values.
#' @param tms list of matrices whose rows are the triple (t,m,sigma) for each band.
#' @param beta vector of intercepts
#' @param a amplitude estimates
#' @param atilde prior vector
#' @param rho_anchor vector of the current estimates of the phase
#' @param omega frequency
#' @param gamma1 nonnegative regularization parameter for shrinking amplitudes
#' @param gamma2 nonnegative regularization parameter for shrinking phases
#' @export
mm_phase_obj <- function(rho,tms,beta,a,atilde,rho_anchor,omega,gamma1,gamma2) {
  B <- length(tms)
  L <- update_Lipschitz(tms,beta,a)
  g <- 0
  for (b in 1:B) {
    band <- tms[[b]]
    n <- nrow(band)
    w <- 1/(band[,3]**2)
    theta <- omega*band[,1] + rho_anchor[b]
    s <- sin(theta)
    c <- cos(theta)
    m <- band[,2]
    r <- m - beta[b] - a[b]*s
    f <- (0.5)*t(r)%*%(w*r)
    df <- a[b]*(t(beta[b] + a[b]*s-m)%*%(w*c))
    g <- g + f + df*(rho[b]-rho_anchor[b]) + 0.5*L[b]*(rho[b]-rho_anchor[b])**2
  }
  #  penalty <- 0
  #  for (b in 1:(B-1)) {
  #    for (bb in (b+1):B) {
  #      penalty <- penalty + 0.5*sum((rho[b] - rho[bb])**2)
  #    }
  #  }
  penalty <- gamma2*(sum(rho**2) - (1/B)*(sum(rho)**2))
  penalty <- penalty + gamma1*(sum(a**2) - (sum(atilde*a))**2)
  return(g + 0.5*penalty)
}

#' Update Phase parameter
#' 
#' \code{update_rho_inexact} inexactly updates the phase parameter rho via an MM algorithm using a convex quadratic majorization.
#' 
#' @param tms list of matrices whose rows are the triple (t,mu,sigma) for each band.
#' @param beta vector of the current intercept estimates
#' @param a amplitude estimates
#' @param rho vector of the current estimates of the phase
#' @param omega frequency
#' @param gamma nonnegative regularization parameter
#' @param max_iter maximum number of iterations
#' @export
#' @examples
#' test_data <- synthetic_multiband()
#' tms <- test_data$tms
#' B <- test_data$B
#' beta <- test_data$beta
#' a <- test_data$a
#' rho <- test_data$rho
#' omega <- test_data$omega
#' gamma <- 1
#' 
#' ## Check answer
#' rho_next <- update_rho_inexact(tms,beta,a,rho,omega,gamma,max_iter=1)
#' 
#' L <- update_Lipschitz(tms,beta,a)
#' f <- L + gamma
#' zeta <- update_zeta(tms,beta,a,rho,L,omega)
#' rho_direct <- solve(diag(f)-(gamma/B),zeta)
#' norm(as.matrix(rho_direct-rho_next),'f')
#' 
#' ## Verify monotonicity of MM algorithm
#' max_iter <- 1e2
#' obj <- double(max_iter)
#' loss <- double(max_iter)
#' rho_last <- rho
#' at <- rep(1/sqrt(B),B)
#' for (iter in 1:max_iter) {
#'   rho_next <- update_rho_inexact(tms,beta,a,rho_last,omega,gamma,max_iter=1)
#'   obj[iter] <- mm_phase_obj(rho_next,tms,beta,a,at,rho_last,omega,gamma,gamma)
#'   loss[iter] <- pnll(tms,beta,a,at,rho_next,omega,gamma,gamma) 
#'   rho_last <- rho_next
#' }
#' obj <- c(mm_phase_obj(rho,tms,beta,a,at,rho,omega,gamma,gamma),obj)
#' plot(1:(max_iter+1),obj,xlab='iteration',ylab='mm objective',pch=16)
#' loss <- c(pnll(tms,beta,a,at,rho,omega,gamma,gamma),loss)
#' plot(1:(max_iter+1),loss,xlab='iteration',ylab='loss',pch=16)
update_rho_inexact <- function(tms,beta,a,rho,omega,gamma,max_iter=5) {
  B <- length(tms)
  L <- update_Lipschitz(tms,beta,a)
  f <- L + gamma
  for (iter in 1:max_iter) {
    zeta <- update_zeta(tms,beta,a,rho,L,omega)
    fizeta <- zeta/f
    rho <- fizeta + gamma*(sum(fizeta)/(B - gamma*sum(1/f)))/f
    ## Wrap the rho step back into the interval [-pi,pi]
    #    print(rho_next)
  }
  rho <- ((rho + pi) %% (2*pi)) - pi
  return(rho)
}
