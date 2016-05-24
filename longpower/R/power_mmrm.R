#' Linear mixed model sample size calculations.
#' 
#' This function performs the sample size calculation for a mixed model of
#' repeated measures with general correlation structure. See Lu, Luo, & Chen
#' (2008) for parameter definitions and other details. This function executes
#' Formula (3) on page 4.
#' 
#' See Lu, Luo, & Chen (2008).
#' 
#' @param N total sample size
#' @param Ra correlation matrix for group a
#' @param ra retention in group a
#' @param sigmaa standard deviation of observation of interest in group a
#' @param Rb correlation matrix for group a
#' @param rb retention in group b
#' @param sigmab standard deviation of observation of interest in group b. If
#' NULL, \code{sigmab} is assumed same as \code{sigmaa}. If not NULL,
#' \code{sigmaa} and \code{sigmab} are averaged.
#' @param lambda allocation ratio
#' @param delta effect size
#' @param sig.level type one error
#' @param power power
#' @param alternative one- or two-sided test
#' @return The number of subject required per arm to attain the specified
#' \code{power} given \code{sig.level} and the other parameter estimates.
#' @author Michael C. Donohue
#' @seealso \code{\link{power.mmrm.ar1}}, \code{\link{lmmpower}},
#' \code{\link{diggle.linear.power}}
#' @references Lu, K., Luo, X., Chen, P.-Y. (2008) Sample size estimation for
#' repeated measures analysis in randomized clinical trials with missing data.
#' \emph{International Journal of Biostatistics}, 4, (1)
#' @keywords power sample size mixed effects random effects
#' @examples
#' 
#' # reproduce Table 1 from Lu, Luo, & Chen (2008)
#' phi1 <- c(rep(1, 6), 2, 2)
#' phi2 <- c(1, 1, rep(2, 6))
#' lambda <- c(1, 2, sqrt(1/2), 1/2, 1, 2, 1, 2)
#' ztest <- ttest1 <- c()
#' for(i in 1:8){
#'   Na <- (phi1[i] + lambda[i] * phi2[i])*(qnorm(0.05/2) + qnorm(1-0.90))^2*(0.5^-2)
#'   Nb <- Na/lambda[i]
#'   ztest <- c(ztest, Na + Nb)
#'   v <- Na + Nb - 2
#'   Na <- (phi1[i] + lambda[i] * phi2[i])*(qt(0.05/2, df = v) + qt(1-0.90, df = v))^2*(0.5^-2)
#'   Nb <- Na/lambda[i]
#'   ttest1 <- c(ttest1, Na + Nb)
#' }
#' data.frame(phi1, phi2, lambda, ztest, ttest1)
#' 
#' Ra <- matrix(0.25, nrow = 4, ncol = 4)
#' diag(Ra) <- 1
#' ra <- c(1, 0.90, 0.80, 0.70)
#' sigmaa <- 1
#' power.mmrm(Ra = Ra, ra = ra, sigmaa = sigmaa, delta = 0.5, power = 0.80)
#' power.mmrm(N = 174, Ra = Ra, ra = ra, sigmaa = sigmaa, delta = 0.5)
#' power.mmrm(N = 174, Ra = Ra, ra = ra, sigmaa = sigmaa, power = 0.80)
#' 
#' power.mmrm(Ra = Ra, ra = ra, sigmaa = sigmaa, delta = 0.5, power = 0.80, lambda = 2)
#' power.mmrm(N = 174, Ra = Ra, ra = ra, sigmaa = sigmaa, delta = 0.5, lambda = 2)
#' power.mmrm(N = 174, Ra = Ra, ra = ra, sigmaa = sigmaa, power = 0.80, lambda = 2)
#' 
#' @export power.mmrm
power.mmrm <- function(N = NULL, Ra = NULL, ra = NULL, sigmaa = NULL, 
  Rb = NULL, rb = NULL, sigmab = NULL, lambda = 1,
  delta = NULL, sig.level = 0.05, power = NULL, 
  alternative = c("two.sided", "one.sided"))
{
  if (sum(sapply(list(N, delta, power, sig.level), is.null)) != 1) 
      stop("exactly one of 'N', 'delta', 'power', and 'sig.level' must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
      sig.level | sig.level > 1)) 
      stop("'sig.level' must be numeric in [0, 1]")
  alternative <- match.arg(alternative)

  # formula (3) on page 4 in
  # Lu, K., Luo, X., & Chen, P.-Y. (July 14, 2008). Sample size estimation for repeated measures analysis in randomized clinical trials with missing data. International Journal of Biostatistics, 4, (1)
    
  n.body <- quote({
    Ia <- 0
    ra <- c(ra, 0)
    for(j in 1:nrow(Ra)){
      Raj <- matrix(0, nrow(Ra), nrow(Ra))
      Raj[1:j, 1:j] <- solve(Ra[1:j, 1:j])
      Ia <- Ia + (ra[j] - ra[j+1]) * Raj
    }
    phia <- solve(Ia)[j,j]

    if(is.null(rb)) rb <- ra
    if(is.null(Rb)) Rb <- Ra
    Ib <- 0
    rb <- c(rb, 0)
    for(j in 1:nrow(Rb)){
      Rbj <- matrix(0, nrow(Rb), nrow(Rb))
      Rbj[1:j, 1:j] <- solve(Rb[1:j, 1:j])
      Ib <- Ib + (rb[j] - rb[j+1]) * Rbj
    }
    phib <- solve(Ib)[j,j]  

    if(is.null(sigmab)){
      sigma <- sigmaa
    }else{
      sigma <- mean(c(sigmaa, sigmab))
    }
  
    Na <- as.numeric(
      (phia + lambda * phib)*(qnorm(ifelse(alternative=="two.sided", sig.level/2, sig.level)) + qnorm(1-power))^2*
      sigma^2/delta^2
    )
    Nb <- as.numeric(Na/lambda)
    list(N = Na + Nb, Na = Na, Nb = Nb)
  })
    
  if (is.null(sig.level)) 
      sig.level <- uniroot(function(sig.level) eval(n.body)$N - 
          N, c(1e-10, 1 - 1e-10))$root
  else if (is.null(power)) 
      power <- uniroot(function(power) eval(n.body)$N - 
          N, c(1e-3, 1 - 1e-10))$root
  else if (is.null(delta)) 
      delta <- uniroot(function(delta) eval(n.body)$N - 
          N, c(1e-10, 1e5))$root
    
  Ns <- eval(n.body)
  Na <- Ns$Na
  Nb <- Ns$Nb
  METHOD <- "Power for Mixed Model of Repeated Measures (Lu, Luo, & Chen, 2008)"
  structure(list(n1 = Na, n2 = Nb, 
        retention1 = ra[-length(ra)], retention2 = rb[-length(rb)], 
        delta = delta, sig.level = sig.level, 
        power = power, alternative = alternative, 
        method = METHOD), class = "power.htest")
}





#' Linear mixed model sample size calculations.
#' 
#' This function performs the sample size calculation for a mixed model of
#' repeated measures with AR(1) correlation structure. See Lu, Luo, & Chen
#' (2008) for parameter definitions and other details.
#' 
#' See Lu, Luo, & Chen (2008).
#' 
#' @param N total sample size
#' @param rho AR(1) correlation parameter
#' @param ra retention in group a
#' @param sigmaa standard deviation of observation of interest in group a
#' @param rb retention in group a (assumed same as \code{ra} if left blank)
#' @param sigmab standard deviation of observation of interest in group b. If
#' NULL, \code{sigmab} is assumed same as \code{sigmaa}. If not NULL,
#' \code{sigmaa} and \code{sigmab} are averaged.
#' @param lambda allocation ratio
#' @param times observation times
#' @param delta effect size
#' @param sig.level type one error
#' @param power power
#' @param alternative one- or two-sided test
#' @return The number of subject required per arm to attain the specified
#' \code{power} given \code{sig.level} and the other parameter estimates.
#' @author Michael C. Donohue
#' @seealso \code{\link{power.mmrm}}, \code{\link{lmmpower}},
#' \code{\link{diggle.linear.power}}
#' @references Lu, K., Luo, X., Chen, P.-Y. (2008) Sample size estimation for
#' repeated measures analysis in randomized clinical trials with missing data.
#' \emph{International Journal of Biostatistics}, 4, (1)
#' @keywords power sample size mixed effects random effects
#' @examples
#' 
#' # reproduce Table 2 from Lu, Luo, & Chen (2008)
#' tab <- c()
#' for(J in c(2,4))
#' for(aJ in (1:4)/10)
#' for(p1J in c(0, c(1, 3, 5, 7, 9)/10)){
#'   rJ <- 1-aJ
#'   r <- seq(1, rJ, length = J)
#'   # p1J = p^(J-1)
#'   tab <- c(tab, power.mmrm.ar1(rho = p1J^(1/(J-1)), ra = r, sigmaa = 1, 
#'     lambda = 1, times = 1:J,
#'     delta = 1, sig.level = 0.05, power = 0.80)$phi1)
#' }
#' matrix(tab, ncol = 6, byrow = TRUE)
#' 
#' # approximate simulation results from Table 5 from Lu, Luo, & Chen (2008)
#' ra <- c(100, 76, 63, 52)/100
#' rb <- c(100, 87, 81, 78)/100
#' 
#' power.mmrm.ar1(rho=0.6, ra=ra, sigmaa=1, rb = rb, 
#'   lambda = sqrt(1.25/1.75), power = 0.904, delta = 0.9
#' )
#' power.mmrm.ar1(rho=0.6, ra=ra, sigmaa=1, rb = rb, 
#'   lambda = 1.25/1.75, power = 0.910, delta = 0.9
#' )
#' power.mmrm.ar1(rho=0.6, ra=ra, sigmaa=1, rb = rb, 
#'   lambda = 1, power = 0.903, delta = 0.9
#' )
#' power.mmrm.ar1(rho=0.6, ra=ra, sigmaa=1, rb = rb, 
#'   lambda = 2, power = 0.904, delta = 0.9
#' )
#' 
#' power.mmrm.ar1(N=81, ra=ra, sigmaa=1, rb = rb, 
#'   lambda = sqrt(1.25/1.75), power = 0.904, delta = 0.9
#' )
#' power.mmrm.ar1(N=87, rho=0.6, ra=ra, sigmaa=1, rb = rb, 
#'   lambda = 1.25/1.75, power = 0.910
#' )
#' power.mmrm.ar1(N=80, rho=0.6, ra=ra, sigmaa=1, rb = rb, 
#'   lambda = 1, delta = 0.9
#' )
#' power.mmrm.ar1(N=84, rho=0.6, ra=ra, sigmaa=1, rb = rb, 
#'   lambda = 2, power = 0.904, delta = 0.9, sig.level = NULL
#' )
#' 
#' @export power.mmrm.ar1
power.mmrm.ar1 <- function(N = NULL, rho = NULL, 
  ra = NULL, sigmaa = NULL, rb = NULL, sigmab = NULL, 
  lambda = 1, times = 1:length(ra),
  delta = NULL, sig.level = 0.05, power = NULL, 
  alternative = c("two.sided", "one.sided"))
{
  if (sum(sapply(list(N, rho, delta, power, sig.level), is.null)) != 1) 
      stop("exactly one of 'N', 'rho', 'delta', 'power', and 'sig.level' must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
      sig.level | sig.level > 1)) 
      stop("'sig.level' must be numeric in [0, 1]")
  alternative <- match.arg(alternative)

  # formula (3) on page 4 in
  # Lu, K., Luo, X., & Chen, P.-Y. (July 14, 2008). Sample size estimation for repeated measures analysis in randomized clinical trials with missing data. International Journal of Biostatistics, 4, (1)
  if(length(times) != length(ra)) stop("ra and times should be the same length")
  
  phia <- phib <- NULL
  
  n.body <- quote({
      J <- length(ra)
      phia <- 1/ra[J] - sum(
        rho^(2*(times[J] - times[-J])) *
        (1/ra[-1] - 1/ra[-J]))

      if(!is.null(rb)){
        phib <- 1/rb[J] - sum(
          rho^(2*(times[J] - times[-J])) *
          (1/rb[-1] - 1/rb[-J]))
      }else{
        rb <- ra
        phib <- phia
      }
      if(is.null(sigmab)){
        sigma <- sigmaa
      }else{
        sigma <- mean(c(sigmaa, sigmab))
      }
  
      Na <- as.numeric(
        (phia + lambda * phib)*(qnorm(ifelse(alternative=="two.sided", sig.level/2, sig.level)) + qnorm(1-power))^2*
        sigma^2/delta^2
      )
      Nb <- as.numeric(Na/lambda)
      list(N = Na + Nb, Na = Na, Nb = Nb)
  })
  
  if (is.null(sig.level)) 
      sig.level <- uniroot(function(sig.level) eval(n.body)$N - 
          N, c(1e-10, 1 - 1e-10))$root
  else if (is.null(power)) 
      power <- uniroot(function(power) eval(n.body)$N - 
          N, c(1e-3, 1 - 1e-10))$root
  else if (is.null(delta)) 
      delta <- uniroot(function(delta) eval(n.body)$N - 
          N, c(1e-10, 1e5))$root
  else if (is.null(rho)) 
      rho <- uniroot(function(rho) eval(n.body)$N - 
          N, c(1e-10, 1 - 1e-10))$root
    
  Ns <- eval(n.body)
  Na <- Ns$Na
  Nb <- Ns$Nb
  
  METHOD <- "Power for Mixed Model of Repeated Measures (Lu, Luo, & Chen, 2008)"
  structure(list(n1 = Na, n2 = Nb, rho = rho, 
        retention1 = ra[-length(ra)], retention2 = rb[-length(rb)],
        phi1 = phia, phi2 = phib, 
        delta = delta, times = times, sig.level = sig.level, 
        power = power, alternative = alternative, 
        method = METHOD), class = "power.htest")
}
