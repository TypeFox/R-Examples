DressEnsemble <- function(ens, dressing.method="silverman", 
                          parameters=NA) 
{

  # preprocess
  l <- Preprocess(ens=ens)
  ens <- l[["ens"]]

  # silverman's rule of thumb
  if (dressing.method == "silverman") {
    n.members <- rowSums(!is.na(ens))
    if (any(n.members==1)) {
      warning(c("Some ensembles have only one member. ",
                "The kernel width is set to zero for these."))
    }
    stdevs <- apply(ens, 1, sd, na.rm=TRUE)
    ker.wd <- (4 * stdevs^5 / (3 * n.members))^0.2

    K <- max(n.members)
    ker.wd <- matrix(rep(ker.wd, K), ncol=K)

    ker.type <- "gauss"
  }

  #
  # affine kernel dressing
  #
  #       p(y|x) = 1 / K * sum {dnorm(y, z.i(x), s(x))}
  # where   s(x) = (4/3/K)^0.4 * (s1 + s2 * a^2 * var(x))
  # and   z.i(x) = r1 + r2 * mean(x) + a * x[i]
  #
  # dressing parameters a, r1, r2, s1, s2 are provided by the user
  if (dressing.method=="akd") {

    stopifnot(all(names(parameters) %in% c("a", "r1", "r2", "s1", "s2")))

    K <- max(rowSums(!is.na(ens)))
    v.x <- apply(ens, 1, var, na.rm=TRUE)
    m.x <- rowMeans(ens, na.rm=TRUE)
    sf.2 <- (4 / 3 / K) ^ 0.4

    z <- with(parameters, r1 + r2 * m.x  + a * ens)
    ens <- z

    s2 <- with(parameters, sf.2 * (s1 + s2 * a * a * v.x))
    s <- sapply(s2, function(x) sqrt(max(x, 0)))
    ker.wd <- matrix(rep(s, K), ncol=K)

    ker.type <- "gauss"
  }


  #
  # affine kernel dressing as above, but parameter are estimated by optimizing
  # crps with respect to obs (which must be provided as parameters)
  #
  if (dressing.method=="akd.fit") {

    stopifnot(any(names(parameters) == "obs"))
    obs <- parameters[["obs"]]

    stopifnot(length(obs) == nrow(ens))

    parms <- as.list(FitAkdParameters(ens, obs))

    K <- max(rowSums(!is.na(ens)))
    v.x <- apply(ens, 1, var, na.rm=TRUE)
    m.x <- rowMeans(ens, na.rm=TRUE)
    sf.2 <- (4 / 3 / K) ^ 0.4

    z <- with(parms, r1 + r2 * m.x  + a * ens)
    ens <- z

    s2 <- with(parms, sf.2 * (s1 + s2 * a * a * v.x))
    s <- sapply(s2, function(x) sqrt(max(x, 0)))
    ker.wd <- matrix(rep(s, K), ncol=K)

    ker.type <- "gauss"
  }

  # create object
  dressed.ens <- list(ens=ens, ker.wd=ker.wd, ker.type=ker.type)
  class(dressed.ens) <- "dressed.ens"

  # return
  dressed.ens
}


