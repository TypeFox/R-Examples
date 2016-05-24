#' @rdname U-utils
#' @export
qU <- function(p, beta, distname, use.mean.variance = TRUE) {

  check_distname(distname)
  names(beta) <- get_beta_names(distname)
  check_beta(beta, distname = distname)
  
  sigma.x <- beta2tau(beta, distname, 
                      use.mean.variance = use.mean.variance)["sigma_x"]
  switch(distname,
         cauchy = {
           qU <- function(p) qcauchy(p)
         },
         chisq = {
           qU <- function(p) qchisq(p, df = beta) / sigma.x
         },
         exp = {
           qU <- function(p) qexp(p)
         },
         "f" = {
           qU <- function(p) qf(p, beta[1], beta[2]) / sigma.x
         },
         gamma = {
           qU <- function(p) qgamma(p, shape = beta["shape"], 
                                    scale = beta["scale"]) / sigma.x
         },
         laplace = {
           #TODO
         },
         normal = {
           qU <- function(p) qnorm(p)
         },
         t = {
           if (beta["df"] <= 2 && use.mean.variance) {
             stop("'df' of t-distribution for location-scale Lambert W x t distributions must be ",
                  " larger than 2.")
           }
           ss <- sigma.x / beta["scale"]
           names(ss) <- NULL
           qU <- function(p) qt(p, df = beta["df"]) / ss
         }, 
         unif = {
           qU <- function(p) qunif(p, -sqrt(3), sqrt(3))
         })
  
  q.u <- qU(p)
  names(q.u) <- names(p)
  return(q.u)
} 