#' @rdname U-utils
#' @export
rU <- function(n, beta, distname, use.mean.variance = TRUE) {
  check_distname(distname)
  names(beta) <- get_beta_names(distname)
  
  sigma.x <- beta2tau(beta, distname = distname,
                      use.mean.variance = use.mean.variance)["sigma_x"]
  switch(distname,
         cauchy = {
           uu <- rcauchy(n)
         },
         chisq = {
           uu <- rchisq(n, df = beta) / sigma.x
         },
         exp = {
           uu <- rexp(n)
         },
         #f = {
        #   #TODO
        # },
         gamma = {
           uu <- rgamma(n, shape = beta["shape"], scale = beta["scale"]) / sigma.x
         },
         laplace = {
           #TODO
         },
         normal = {
           uu <- rnorm(n)
         },
         t = {
           ss<- sigma.x / beta["scale"]
           uu <- rt(n, df = beta["df"]) / ss
         }, 
         unif = {
           uu <- runif(n, -sqrt(3), sqrt(3))
         }
  )
  return(uu)
} 
