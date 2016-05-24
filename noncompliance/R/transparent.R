#' @useDynLib noncompliance
NULL

###############################################################################
# Bayesian procedures for transparent reparametrization
###############################################################################

#' Check of the Instrumental Variable (IV) inequalities.
#'
#' This checks whether the Instrumental Variable (IV) inequalities
#'    for a binary dataset have been \emph{satisfied} empirically,
#'    assuming only Randomization and Exclusion Restriction for the
#'    principal strata of Always Takers and Never Takers.
#'    Monotonicity (no Defiers) is \emph{not} assumed here.
#'
#' @export
#' @param n_y0x0z0 Number of individuals with Y=0, X=0, Z=0.
#'    Alternatively, a vector with elements
#'    (either counts, p(y, x , z) or p(y, x | z)) in the order of the arguments.
#' @param n_y1x0z0 Number of individuals with Y=1, X=0, Z=0.
#' @param n_y0x1z0 Number of individuals with Y=0, X=1, Z=0.
#' @param n_y1x1z0 Number of individuals with Y=1, X=1, Z=0.
#' @param n_y0x0z1 Number of individuals with Y=0, X=0, Z=1.
#' @param n_y1x0z1 Number of individuals with Y=1, X=0, Z=1.
#' @param n_y0x1z1 Number of individuals with Y=0, X=1, Z=1.
#' @param n_y1x1z1 Number of individuals with Y=1, X=1, Z=1.
#' @param verbose Whether to return all the IV inequalities (TRUE)
#'    or just a check that the inequalities have been satisfied empirically.
#' @return A list of all the IV inequalities or a check of whether all
#'    the inequalities have been satisfied empirically.
#' @examples
#' Check_IV_ineqs(158, 14, 0, 0, 52, 12, 23, 78)
#' Check_IV_ineqs(c(158, 14, 0, 0, 52, 12, 23, 78))
#' Check_IV_ineqs(158, 14, 0, 0, 52, 12, 23, 78, TRUE)
#' Check_IV_ineqs(99, 1027, 30, 233, 84, 935, 31, 422)
#' Check_IV_ineqs(c(99, 1027, 30, 233, 84, 935, 31, 422))
#' Check_IV_ineqs(99, 1027, 30, 233, 84, 935, 31, 422, TRUE)
#' @references {A. Balke and J. Pearl. (1997).
#'    Bounds on treatment effects from studies with imperfect compliance.
#'    \emph{Journal of the American Statistical Association, 1171-1176}.},
#' @references {B. Bonet. (2001).
#'    Instrumentality tests revisited.
#'    \emph{In Proceedings of the Seventeenth Conference on Uncertainty in
#'    Artificial Intelligence, 48-55}.}

Check_IV_ineqs <- function(n_y0x0z0, n_y1x0z0=NA, n_y0x1z0=NA, n_y1x1z0=NA,
                           n_y0x0z1=NA, n_y1x0z1=NA, n_y0x1z1=NA, n_y1x1z1=NA,
                           verbose=FALSE){

  if (length(n_y0x0z0) == 8) {
    n_yxz <- n_y0x0z0; n_y0x0z0 <- NULL
    n_y0x0z0 <- n_yxz[1]; n_y1x0z0 <- n_yxz[2]
    n_y0x1z0 <- n_yxz[3]; n_y1x1z0 <- n_yxz[4]
    n_y0x0z1 <- n_yxz[5]; n_y1x0z1 <- n_yxz[6]
    n_y0x1z1 <- n_yxz[7]; n_y1x1z1 <- n_yxz[8]
  }
  pyx.z0 <- c(n_y0x0z0, n_y1x0z0, n_y0x1z0, n_y1x1z0)
  pyx.z0 <- pyx.z0 / sum(pyx.z0)
  pyx.z1 <- c(n_y0x0z1, n_y1x0z1, n_y0x1z1, n_y1x1z1)
  pyx.z1 <- pyx.z1 / sum(pyx.z1)

  ivineq <- c( pyx.z0[1] + pyx.z1[2], pyx.z0[2] + pyx.z1[1],  #4
               pyx.z0[3] + pyx.z1[4], pyx.z0[4] + pyx.z1[3] ) #5

  if( any( ivineq > 1 ) ){
    violated <- "p(y, x | z) does not satisfy inequalities for IV Model!"
  } else{
    violated <- "All IV inequalities satisfied: distribution compatible with IV model."
  }

  if(verbose){
    ivineq <- round(ivineq, 4)
    return( list( list(
      "p(y0, x0 | z0) + p(y1, x0 | z1)" = ivineq[1],
      "p(y1, x0 | z0) + p(y0, x0 | z1)" = ivineq[2],
      "p(y0, x1 | z0) + p(y1, x1 | z1)" = ivineq[3],
      "p(y1, x1 | z0) + p(y0, x1 | z1)" = ivineq[4]),
      violated ))
  } else {
    return( all( ivineq <= 1 ) )
  }
}

#' Bounds for the Average Causal Effect (ACE).
#'
#' The empirical bounds for the Average Causal Effect (ACE),
#'    under the assumptions of the Instrumental Variable (IV) model.
#'
#' @export
#' @inheritParams Check_IV_ineqs
#' @return The empirical bounds for the ACE.
#' @examples
#' ACE_bounds(158, 14, 0, 0, 52, 12, 23, 78)
#' ACE_bounds(c(158, 14, 0, 0, 52, 12, 23, 78))
#' ACE_bounds(99, 1027, 30, 233, 84, 935, 31, 422)
#' ACE_bounds(c(99, 1027, 30, 233, 84, 935, 31, 422))
#' @references {Richardson, T. S.; Robins, J. M. (2014).
#'    ACE Bounds; SEMs with Equilibrium Conditions.
#'    \emph{Statist. Sci. 29, no. 3, 363-366.}.}

ACE_bounds <- function(n_y0x0z0, n_y1x0z0=NA, n_y0x1z0=NA, n_y1x1z0=NA,
                       n_y0x0z1=NA, n_y1x0z1=NA, n_y0x1z1=NA, n_y1x1z1=NA){

  if (length(n_y0x0z0) == 8) {
    n_yxz <- n_y0x0z0; n_y0x0z0 <- NULL
    n_y0x0z0 <- n_yxz[1]; n_y1x0z0 <- n_yxz[2]
    n_y0x1z0 <- n_yxz[3]; n_y1x1z0 <- n_yxz[4]
    n_y0x0z1 <- n_yxz[5]; n_y1x0z1 <- n_yxz[6]
    n_y0x1z1 <- n_yxz[7]; n_y1x1z1 <- n_yxz[8]
  }

  pyx.z <- expand.grid(X=0:1, Y=0:1, Z=0:1, p=NA)
  pyx.z[pyx.z$Y == 0 & pyx.z$X == 0 & pyx.z$Z == 0, "p"] <- n_y0x0z0
  pyx.z[pyx.z$Y == 1 & pyx.z$X == 0 & pyx.z$Z == 0, "p"] <- n_y1x0z0
  pyx.z[pyx.z$Y == 0 & pyx.z$X == 1 & pyx.z$Z == 0, "p"] <- n_y0x1z0
  pyx.z[pyx.z$Y == 1 & pyx.z$X == 1 & pyx.z$Z == 0, "p"] <- n_y1x1z0
  pyx.z[pyx.z$Y == 0 & pyx.z$X == 0 & pyx.z$Z == 1, "p"] <- n_y0x0z1
  pyx.z[pyx.z$Y == 1 & pyx.z$X == 0 & pyx.z$Z == 1, "p"] <- n_y1x0z1
  pyx.z[pyx.z$Y == 0 & pyx.z$X == 1 & pyx.z$Z == 1, "p"] <- n_y0x1z1
  pyx.z[pyx.z$Y == 1 & pyx.z$X == 1 & pyx.z$Z == 1, "p"] <- n_y1x1z1

  for (z in 0:1) {
    pyx.z[pyx.z$Z == z, "p"] <- pyx.z[pyx.z$Z == z, "p"] /
      sum( pyx.z[pyx.z$Z == z, "p"] )
  }

  g_ij <- expand.grid(X=0:1, Y=0:1, g=NA)
  for (i in 0:1) {
    for (j in 0:1) {
      t1 <- sapply(0:1,function(z) {
        pyx.z[pyx.z$X == i & pyx.z$Y == j & pyx.z$Z == z, "p"] +
          sum( pyx.z[pyx.z$X == (1-i) & pyx.z$Z == z, "p"] )
      })
      t2 <- sapply(0:1,function(z) {
        pyx.z[pyx.z$X == i & pyx.z$Y == j & pyx.z$Z == z, "p"] +
          pyx.z[pyx.z$X == (1-i) & pyx.z$Y == 0 & pyx.z$Z == z, "p"] +
          pyx.z[pyx.z$X == i & pyx.z$Y == j & pyx.z$Z == (1-z), "p"] +
          pyx.z[pyx.z$X == (1-i) & pyx.z$Y == 1 & pyx.z$Z == (1-z), "p"]
      })
      g_ij[g_ij$X == i & g_ij$Y == j, "g"] <- min( min(t1), min(t2) )
    }
  }

  ace.lb <- 1
  ace.ub <- -1
  for (i in 0:1){
    ace.lb <- ace.lb - g_ij[g_ij$X == i & g_ij$Y == (1-i), "g"]
    ace.ub <- ace.ub + g_ij[g_ij$X == i & g_ij$Y == i, "g"]
  }
  return(c(ace.lb, ace.ub))
}

#' Posterior bounds for the Average Causal Effect (ACE).
#'
#' The posterior bounds for the Average Causal Effect (ACE) is found
#'    based on a transparent reparametrization (see reference below),
#'    using a Dirichlet prior.
#'    A binary Instrumental Variable (IV) model is assumed here.
#'
#' @export
#' @importFrom stats rgamma
#' @inheritParams Check_IV_ineqs
#' @param n_y0x0z0 Number of individuals with Y=0, X=0, Z=0.
#'    Alternatively, a vector with elements in the order of the arguments.
#' @param prior Hyperparameters for the Dirichlet prior for p(y, x | z),
#'    in the order of the arguments.
#' @param num.sims Number of Monte Carlo draws from the posterior.
#' @return A data frame with the posterior bounds for the ACE,
#'    based only on sampled distributions (from the posterior)
#'    that satisfied the IV inequalites.
#' @examples
#' ACE_bounds_posterior(158, 14, 0, 0, 52, 12, 23, 78,
#'      prior = c( rep(1, 2), rep(0, 2), rep(1, 4)))
#' ACE_bounds_posterior(158, 14, 0, 0, 52, 12, 23, 78,
#'      prior = c( rep(1/2, 2), rep(0, 2), rep(1/4, 4)))
#' \dontrun{
#' ace.bnds.lipid <- ACE_bounds_posterior(158, 14, 0, 0, 52, 12, 23, 78,
#'      prior = c( rep(1, 2), rep(0, 2), rep(1, 4)), num.sims = 2e4)
#' summary(ace.bnds.lipid) }
#' @references {Richardson, T. S., Evans, R. J., & Robins, J. M. (2011).
#'    Transparent parameterizations of models for potential outcomes.
#'    \emph{Bayesian Statistics, 9, 569-610}.}

ACE_bounds_posterior <- function(
  n_y0x0z0, n_y1x0z0=NA, n_y0x1z0=NA, n_y1x1z0=NA,
  n_y0x0z1=NA, n_y1x0z1=NA, n_y0x1z1=NA, n_y1x1z1=NA,
  prior, num.sims = 1e3 ) {

  if (length(n_y0x0z0) == 8) {
    n_yxz <- n_y0x0z0; n_y0x0z0 <- NULL
    n_y0x0z0 <- n_yxz[1]; n_y1x0z0 <- n_yxz[2]
    n_y0x1z0 <- n_yxz[3]; n_y1x1z0 <- n_yxz[4]
    n_y0x0z1 <- n_yxz[5]; n_y1x0z1 <- n_yxz[6]
    n_y0x1z1 <- n_yxz[7]; n_y1x1z1 <- n_yxz[8]
  }
  nyx.z0 <- c(n_y0x0z0, n_y1x0z0, n_y0x1z0, n_y1x1z0)
  nyx.z1 <- c(n_y0x0z1, n_y1x0z1, n_y0x1z1, n_y1x1z1)
  if (sum(nyx.z0) < 1 || sum(nyx.z1) < 1) {
    print("Please check the counts entered.")
  }

  # Posterior
  post.z0 <- prior[1:4] + nyx.z0
  post.z1 <- prior[5:8] + nyx.z1

  # Single posterior draw
  dirichlet <- function(alpha){
    theta <- sapply(alpha, function(a) {
      out <- 0.0
      if(a>0) {
        out <- rgamma(n = 1, shape = a, rate = 1)
      }
      return(out)
    })
    return( theta / sum(theta) )
  }

  # MC samples from posterior
  theta.sims.z0 <- t(replicate( num.sims, dirichlet(post.z0) ))
  theta.sims.z1 <- t(replicate( num.sims, dirichlet(post.z1) ))
  theta.sims <- cbind(theta.sims.z0, theta.sims.z1)

  # Remove sampled distributions that violate the IV inequalites
  all.ivs.ok <- apply(theta.sims, 1, Check_IV_ineqs)
  posterior.theta.sims.iv <- theta.sims[all.ivs.ok, ]

  # Posterior ACE bounds
  ace.bnds.lipid <- t(apply(posterior.theta.sims.iv, 1, ACE_bounds))

  return(ace.bnds.lipid)

}


#' "Triangle" plot of the posterior bounds for the Average Causal Effect (ACE).
#'
#' Plot of the posterior upper bound for the Average Causal Effect (ACE)
#'    against the corresponding lower bound.
#'
#' @export
#' @import graphics
#' @param bounds Posterior bounds from the ACE_bounds_posterior function.
#' @param title.txt Title for the plot.
#' @return A "triangle" plot.
#' @examples
#' ace.bnds.lipid <- ACE_bounds_posterior(158, 14, 0, 0, 52, 12, 23, 78,
#'      prior = c( rep(1, 2), rep(0, 2), rep(1, 4)))
#' ACE_bounds_triangle.plot(ace.bnds.lipid, "Bounds on ACE for Lipid Data")
#' \dontrun{
#' ace.bnds.lipid <- ACE_bounds_posterior(158, 14, 0, 0, 52, 12, 23, 78,
#'      prior = c( rep(1, 2), rep(0, 2), rep(1, 4)), num.sims = 2e4)
#' ACE_bounds_triangle.plot(ace.bnds.lipid, "Bounds on ACE for Lipid Data") }

ACE_bounds_triangle.plot <- function(bounds, title.txt){
  plot( c(-1.001,1.001), c(-1.001,1.001), axes=FALSE, type="n",
        xlab=expression(paste("Lower Bound on ACE")),
        ylab=expression(paste("Upper Bound on ACE")),
        asp=1 )
  axis(1,at=seq(-1,1,by=0.2),cex.axis=1)
  axis(2,at=seq(-1,1,by=0.2),cex.axis=1)
  title(title.txt,cex=3,line=1)
  points(c(-1,1,-1,-1),c(-1,1,1,-1), type="l",lty=1)
  points(bounds[,1], bounds[,2], pch=19, col="#00000110",
         # scale the size of the points based on the number of points
         cex=3/log10(nrow(bounds)))
  points(c(-1,0,0),c(0,0,1),col="red",lty=1,type="l")
}



#' Bounds for the Average Controlled Direct Effect (ACDE).
#'
#' The empirical bounds for the Average Controlled Direct Effect (ACDE) within
#'    the principal strata of Always Takers and Never Takers,
#'    under the assumption of monotonicity (no Defiers).
#'    These are equivalent to an empirical check of the Instrumental Variable
#'    (IV) inequalities (see references below).
#'
#' @export
#' @inheritParams Check_IV_ineqs
#' @param iv.ineqs Whether to return the empirical bounds or
#'    the IV inequalities (TRUE).
#' @return The empirical bounds for the ACDE
#'    among Always Takers and Never Takers, or the
#'    empirical IV inequalities.
#' @examples
#' Check_ACDE_bounds(99, 1027, 30, 233, 84, 935, 31, 422)
#' Check_ACDE_bounds(c(99, 1027, 30, 233, 84, 935, 31, 422))
#' Check_ACDE_bounds(99, 1027, 30, 233, 84, 935, 31, 422, iv.ineqs=TRUE)
#' Check_ACDE_bounds(c(99, 1027, 30, 233, 84, 935, 31, 422), iv.ineqs=TRUE)
#' @references {Richardson, T. S., Evans, R. J., & Robins, J. M. (2011).
#'    Transparent parameterizations of models for potential outcomes.
#'    \emph{Bayesian Statistics, 9, 569-610}.}
#' @references {A. Balke and J. Pearl. (1997).
#'    Bounds on treatment effects from studies with imperfect compliance.
#'    \emph{Journal of the American Statistical Association, 1171-1176}.}

Check_ACDE_bounds <- function(n_y0x0z0, n_y1x0z0=NA, n_y0x1z0=NA, n_y1x1z0=NA,
                              n_y0x0z1=NA, n_y1x0z1=NA, n_y0x1z1=NA, n_y1x1z1=NA,
                              iv.ineqs=FALSE){

  if (length(n_y0x0z0) == 8) {
    n_yxz <- n_y0x0z0; n_y0x0z0 <- NULL
    n_y0x0z0 <- n_yxz[1]; n_y1x0z0 <- n_yxz[2]
    n_y0x1z0 <- n_yxz[3]; n_y1x1z0 <- n_yxz[4]
    n_y0x0z1 <- n_yxz[5]; n_y1x0z1 <- n_yxz[6]
    n_y0x1z1 <- n_yxz[7]; n_y1x1z1 <- n_yxz[8]
  }

  nz0 <- n_y0x0z0+n_y1x0z0+n_y0x1z0+n_y1x1z0
  nz1 <- n_y0x0z1+n_y1x0z1+n_y0x1z1+n_y1x1z1
  nx0.z1 <- max(n_y0x0z1+n_y1x0z1, .Machine$double.eps)
  nx1.z0 <- max(n_y0x1z0+n_y1x1z0, .Machine$double.eps)

  if (iv.ineqs==FALSE) {

    return( list(
      ACDE_NT0.lb = n_y1x0z1/nx0.z1 - min(1, (n_y1x0z0/nz0)/(nx0.z1/nz1)),
      ACDE_NT0.ub = n_y1x0z1/nx0.z1 - max(0, 1 - (n_y0x0z0/nz0)/(nx0.z1/nz1)),
      ACDE_AT1.lb = max(0, 1 - (n_y0x1z1/nz1)/(nx1.z0/nz0)) - n_y1x1z0/nx1.z0,
      ACDE_AT1.ub = min(1, (n_y1x1z1/nz1)/(nx1.z0/nz0)) - n_y1x1z0/nx1.z0 ))

  } else {

    return( list(
      "p(y0,x0|z0) - p(y0,x0|z1)" = (n_y0x0z0/nz0) - (n_y0x0z1/nz1), #12
      "p(y1,x0|z0) - p(y1,x0|z1)" = (n_y1x0z0/nz0) - (n_y1x0z1/nz1), #12
      "p(y1,x1|z1) - p(y1,x1|z0)" = (n_y1x1z1/nz1) - (n_y1x1z0/nz0), #13
      "p(y0,x1|z1) - p(y0,x1|z0)" = (n_y0x1z1/nz1) - (n_y0x1z0/nz0) )) #13

  }
}
