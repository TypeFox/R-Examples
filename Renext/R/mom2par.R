
mom2par <- function(densfun = "exponential",
                    mean, sd = NA) {
  
  if (densfun == "exponential") {

    L <- list(rate = 1/mean)

  } else if (densfun == "gpd") {

    CV <- sd / mean
    xi <- (1- 1/CV/CV)/2
    sigma <- (1-xi)*mean
    L <- list(shape = xi, scale = sigma, loc = 0)

  } else if (densfun == "weibull") {
    n <- max(c(length(mean), length(sd)))
    mean <- rep(mean, length.out = n)
    sd <- rep(sd, length.out = n)
    shape <- rep(NA, n)
    scale <- rep(NA, n)
    ## On trouve le parametre de forme a partir du CV
    ## On vectorise la fonction...
    CVmatch <- function(alphai) {
      sdi / meani - sqrt(gamma(1+2/alphai)/(gamma(1+1/alphai)^2) -1)
    }
    ## use sapply ???
    for (i in 1:n) {
      sdi <- sd[i]
      meani <- mean[i]
      resi <- uniroot(f = CVmatch, interval = c(0.05, 30))
      shape[i] = resi$root
      scale[i] = meani / gamma(1+1/resi$root)
    }
    L <- list(shape = shape, scale = scale)
  } else if (densfun == "gamma") {
     CV <- sd/mean
     shape <- 1/CV/CV
     scale <- mean/shape
     L <- list(shape = shape, scale = scale)
  } else if (densfun == "negative binomial") {
    
    prob <- mean/sd/sd
    size <- mean*prob/(1-prob)
    cat("prob = ", prob, "\n")
    if (prob>1) stop("sd doit etre <= sqrt(mean)")
    L <- list(size = size, prob = prob)

  } else stop("distribution non-reconnue")

  L 

}

