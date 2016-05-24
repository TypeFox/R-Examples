#' Compute empirical variance of encounter rate
#'
#' Computes one of a series of possible variance estimates for the observed
#' encounter rate for a set of sample measurements (e.g., line lengths) and
#' number of observations per sample.
#'
#' The choice of type follows the notation of Fewster et al. (2009) in that
#' there are 8 choices of encounter rate variance that can be computed:
#'
#' \describe{
#' \item{\code{R2}}{random line placement with unequal line lengths (design-assisted estimator)}
#' \item{\code{R3}}{random line placement, model-assisted estimator, based on true contagion process}
#' \item{\code{R4}}{random line placement, model-assisted estimator, based on apparent contagion process}
#' \item{\code{S1}}{systematic line placement, post-stratification with no strata overlap}
#' \item{\code{S2}}{systematic line placement, post-stratification with no strata overlap, variances weighted by line length per stratum}
#' \item{\code{O1}}{systematic line placement, post-stratification with overlapping strata (akin to S1)}
#' \item{\code{O2}}{systematic line placement, post-stratification with overlapping strata (weighted by line length per stratum, akin to S2)}
#' \item{\code{O3}}{systematic line placement, post-stratification with overlapping strata, model-assisted estimator with trend in encounter rate with line length}}
#'
#' Default value is R2, shown in Fewster et al. (2009) to have good performance
#' for completely random designs.  For systematic parallel line transect
#' designs, Fewster et al. recommend O2.
#'
#' For the systematic estimators, pairs are assigned in the order they are
#' given in the \code{lengths} and \code{groups} vectors.
#'
#' @usage   varn(lvec,nvec,type)
#'
#'          covn(lvec, groups1, groups2, type)
#' @aliases varn covn
#' @param lvec vector of sample measurements (e.g., line lengths)
#' @param nvec vector of number observed
#' @param groups1 vector of number of groups observed
#' @param groups2 vector of number of individuals observed
#' @param type choice of variance estimator to use for encounter rate
#' @return Variance of encounter rate as defined by arguments
#' @note This function is also used with different calling arguments to compute
#'   Innes et al variance of the estimated abundances/length rather than
#'   observation encounter rate. The function covn is probably only valid for
#'   R3 and R2.  Currently, the R2 form is used for all types other than R3.
#' @author Jeff Laake
#' @references Fewster, R.M., S.T. Buckland, K.P. Burnham, D.L. Borchers, P.E.
#'   Jupp, J.L. Laake and L. Thomas. 2009. Estimating the encounter rate
#'   variance in distance sampling. Biometrics 65: 225-236.
#' @keywords utility
varn <- function(lvec,nvec,type){
  #  Function courtesy of Rachel Fewster with a few minor changes

  ntot <- sum(nvec)
  L <- sum(lvec)
  k <- length(lvec)

  ## Go through the estimators one by one.
  ## R2, R3, R4, S1, S2, O1, O2, O3
  if(!(type %in% c("R2","R3","R4","S1","S2","O1","O2","O3")))
    stop (paste("Encounter rate variance type '",type,"' is not recognized.",sep=""))

  ## First the estimators based on the assumption of a random sample
  ## of lines: R2, R3, R4:

  ## Estimator R2: var.R2, sd.R2
  if(type=="R2"){
    var.R2 <- (k * sum(lvec^2 * (nvec/lvec - ntot/L)^2))/(L^2 * (k -1))
    return(var.R2)
  }

  ## Estimator R3: var.R3, sd.R3
  if(type=="R3"){
    var.R3 <- 1/(L * (k - 1)) * sum(lvec * (nvec/lvec - ntot/L)^2)
    return(var.R3)
  }

  ## Estimator R4: var.R4, sd.R4
  ## Using approximation to Negbin variance (the apparent contagion model)
  ## using phi = 2 + eps
  if(type=="R4"){
    if(all(lvec==mean(lvec))){
      phi <- (2-2/k)/(1-2/k)
    }else{
      S <- sum(lvec^2)
      C <- sum(lvec^3)
      logvec <- log(lvec)
      D1 <- sum(lvec * logvec)
      D2 <- sum(lvec^2 * logvec)
      D3 <- sum(lvec^3 * logvec)
      eps.top <- 2 * (S^2 - L * C)
      eps.bottom <- L * S * D1 - 2 * S * D2 + 2 * L * D3 - L^2 * D2
      eps <- eps.top/eps.bottom
      phi <- 2 + eps
    }
    alpha <- 1/sum(lvec^phi * (L/lvec - 1))
    var.R4 <- alpha * sum(lvec^phi * (nvec/lvec - ntot/L)^2)
    return(var.R4)
  }

  ## Now the stratified estimators with non-overlapping strata:
  ## S1 and S2:

  ## First group the lines into strata, so that all strata have
  ## two lines but if the last stratum has three if necessary:
  H <- floor(k/2)
  k.h <- rep(2, H)
  if(k %% 2 > 0){
    k.h[H] <- 3
  }
  end.strat <- cumsum(k.h)
  begin.strat <- cumsum(k.h) - k.h + 1

  ## Estimators S1 and S2: var.S1, sd.S1 ; var.S2, sd.S2
  if(type=="S1" | type=="S2"){
    sum.S1 <- 0
    sum.S2 <- 0
    for(h in 1:H){
      nvec.strat <- nvec[begin.strat[h]:end.strat[h]]
      lvec.strat <- lvec[begin.strat[h]:end.strat[h]]
      nbar.strat <- mean(nvec.strat)
      lbar.strat <- mean(lvec.strat)

      ## S1 calculations:
      inner.strat.S1 <- sum((nvec.strat - nbar.strat - (ntot/L) *
        (lvec.strat - lbar.strat))^2)
      sum.S1 <- sum.S1 + k.h[h]/(k.h[h] - 1) * inner.strat.S1

      ## S2 calculations: note that we use estimator R2 within
      ## each stratum:
      L.strat <- sum(lvec.strat)
      var.strat.S2 <- k.h[h]/(L.strat^2 * (k.h[h] - 1)) * sum(
        lvec.strat^2 * (nvec.strat/lvec.strat - nbar.strat/lbar.strat)^2)
      sum.S2 <- sum.S2 + L.strat^2 * var.strat.S2
    }
    if(type=="S1"){
      var.S1 <- sum.S1/L^2
      return(var.S1)
    }else{
      var.S2 <- sum.S2/L^2
      return(var.S2)
    }
  }

  ## Now the stratified estimators with overlapping strata:
  ## O1, O2, O3:
  lvec.1 <- lvec[ - k]
  lvec.2 <- lvec[-1]
  nvec.1 <- nvec[ - k]
  nvec.2 <- nvec[-1]
  ervec.1 <- nvec.1/lvec.1
  ervec.2 <- nvec.2/lvec.2

  ## Estimator O1: var.O1, sd.O1
  if(type=="O1"){
    overlap.varterm <- (nvec.1 - nvec.2 - ntot/L * (lvec.1 - lvec.2))^2
    var.O1 <- k/(2 * L^2 * (k - 1)) * sum(overlap.varterm)
    return(var.O1)
  }

  ## Estimator O2: var.O2, sd.O2
  if(type=="O2"){
    V.overlap.R2 <- ((lvec.1 * lvec.2)/(lvec.1 + lvec.2))^2 * (ervec.1-ervec.2)^2
    var.O2 <- (2 * k)/(L^2 * (k - 1)) * sum(V.overlap.R2)
    return(var.O2)
  }

  ## Estimator O3: var.O3, sd.O3
  if(type=="O3"){
    V.overlap.R3 <- ((lvec.1 * lvec.2)/(lvec.1 + lvec.2)) * (ervec.1 -ervec.2)^2
    var.O3 <- 1/(L * (k - 1)) * sum(V.overlap.R3)
    return(var.O3)
  }
}

covn <- function (lvec, groups1, groups2, type){
  # covn - computes covariance of encounter rate of clusters and individuals as
  #        called from dht.  modeled after varn.
  L <- sum(lvec)
  n1 <- sum(groups1)
  er1 <- n1/L
  n2 <- sum(groups2)
  er2 <- n2/L

  k=length(lvec)
  if(type=="R3"){
    varer <- sum(lvec * (groups1/lvec - er1)*(groups2/lvec-er2))/(L*(k-1))
  }else{
    varer <- k*sum(lvec^2 * (groups1/lvec - er1)*(groups2/lvec-er2))/(L^2*(k-1))
  }

  return(varer)
}
