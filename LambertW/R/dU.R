#' @rdname U-utils
#' @export
dU <- function(u, beta, distname, use.mean.variance = TRUE) {
  check_distname(distname)
  names(beta) <- get_beta_names(distname)
  check_beta(beta, distname)
  
  sigma.x <- beta2tau(beta, distname, 
                      use.mean.variance = use.mean.variance)["sigma_x"]
  
  switch(distname,
         cauchy = {
           fU <- function(u) dcauchy(u)
         }, 
         chisq = {
           fU <- function(u) dchisq(u * sigma.x, df = beta) * sigma.x
         },
         exp = {
           fU <- function(u) dexp(u, rate = 1)
         },
         "f" = {
           fU <- function(u) df(u * sigma.x, beta[1], beta[2]) * sigma.x
         },
         gamma = {
           fU <- function(u) dgamma(u * sigma.x, shape = beta["shape"], 
                                    scale = beta["scale"]) * sigma.x
         },
         # laplace = {
         #   fU = function(u) dlaplace(u, 0, 1/sqrt(2))
         # },
         normal = {
           fU <- function(u) dnorm(u)
         },
         t = {
           ss <- sigma.x / beta["scale"]
           names(ss) <- NULL
           fU <- function(u) dt(u * ss, df = beta["df"]) * ss
         },
         unif = {
           fU <- function(u) dunif(u, -sqrt(3), sqrt(3))
         })
  f.u <- fU(u)
  names(f.u) <- names(u)
  return(f.u)
} 