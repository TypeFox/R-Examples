
## functions for fitting spatial and non-spatial mixture models to fMRI datasets
## the mixture model considered is a mixture of a standard normal distribution and two Gamma functions, denoted N2G
## x ~ p1 * N(0, 1) + p2 * Gamma(a, b) + (1 - p1 -p2) * -Gamma(c, d)
## par = c(a, b, c, d, p1, p2)


N2G <- function(data, par.start = c(4, 2, 4, 2, 0.9, 0.05)) {

    ## fits the N2G model to data using par.start as the starting point of the optimisation
    ## lims is the interval in which observations are assigned to the Normal component
    
    data <- data[data != 0]
    fit <- N2G.Fit(data, par.start, maxit = 500, method = "BFGS")
    lims <- N2G.Region(fit)
    
    return(list(par = fit, lims = lims))
}

N2G.Transform <- function(par) {
    
    ## transform parameters of N2G model so as to lie on the real line
    
    q <- par
    q[1:4] <- log(par[1:4])
    q[5] <- log(par[5] / (1 - par[5] - par[6]))
    q[6] <- log(par[6] / (1 - par[5] - par[6]))
    if(length(par) == 7) q[7] <- log(par[7])
    return(q)
}

N2G.Inverse <- function(par) {
    
    ## transform parameters back to their real domains
    
    q <- par
    q[1:4] <- exp(par[1:4])
    q[5] <- exp(par[5]) / (1 + exp(par[5]) + exp(par[6]))
    q[6] <- exp(par[6]) / (1 + exp(par[5]) + exp(par[6]))
    if(length(par) == 7) q[7] <- exp(par[7])
    
    return(q)
}

N2G.Density <- function (data, par) {
    
    ## density function for the N2G model
    
    pos <- data > 0
    neg <- data < 0
    d <- par[5] * dnorm(data)
    if(par[1] > 0) d[pos] <- d[pos] + par[6] * dgamma(data[pos], par[1], par[2])
    if(par[3] > 0) d[neg] <- d[neg] + (1 - par[5] - par[6]) * dgamma(-data[neg], 
        par[3], par[4])
    return(d)
    
}


N2G.Likelihood <- function(inv.par, data) {

    ## (Negative) Likelihood of the N2G model
    par <- N2G.Inverse(inv.par)
    lik <- -sum(log(N2G.Density(data, par)))

    return(lik)
}

N2G.Class.Probability <- function(data, par){

    ## Posterior Probability of data points being in each class 

    res <- matrix(0, length(data), 3)
    pos <- data > 0
    neg <- data < 0
    res[, 1] <- par[5] * dnorm(data)
    res[pos, 2] <- par[6] * dgamma(data[pos], par[1], par[2])
    res[neg, 3] <- (1 - par[5] - par[6]) * dgamma(-data[neg], par[3], par[4])

    rs <- rowSums(res)
    res <- res / rs
    
    return(res)
}

N2G.Fit <- function(data, par.start, maxit, method){

    ## main fitting function for N2G model
    
    inv.par <- N2G.Transform(par.start)

    ans <- optim(inv.par, fn = N2G.Likelihood, method = method, data = data, control = list(maxit = maxit, trace = 0))
    
    if(ans$convergence!=0) return("No convergence")
    res <- N2G.Inverse(ans$par)

    return(res)
}

N2G.Region <- function(par1) {

    ## calculates the interval within which observations are classified as belonging to the Normal component
    
    f1 <- function(x, p) (p[5] * dnorm(x) - p[6] * dgamma(x, p[1], p[2]))^2
    f2 <- function(x, p) (p[5] * dnorm(x) - (1 - p[5] - p[6]) * dgamma(-x, p[3], p[4]))^2

    ans1 <- optim(par = 1, fn = f1, method = "L-BFGS-B", p = par1, lower = 0, upper = 6, control = list(maxit = 500, trace = 0))$par
    ans2 <- optim(par = -1, fn = f2, method = "L-BFGS-B", p = par1, lower = -6, upper = 0, control = list(maxit = 500, trace = 0))$par
    
    return(c(ans1, ans2))
}



N2G.Likelihood.Ratio <- function(data, par) {

    ## calculates the ratio of the likelihood that data came from the positive Gamma distribution (activation) to the likelihood that data came from the other two distributions (Normal and negative Gamma)
    
    d1 <- array(0, dim = dim(data))
    d0 <- array(0, dim = dim(data))

    p.pos <- par[6]
    p.0   <- par[5]
    p.neg <- 1 - par[5] - par[6]
    
    pos <- data > 0
    neg <- data < 0
    
    d1[pos] <- dgamma(data[pos], par[1], par[2])

    d0 <- (p.0 / (1 - p.pos)) * dnorm(data)
    d0[neg] <- d0[neg] + (p.neg / (1 - p.pos)) * dgamma(-data[neg], par[3], par[4])

    ans <- d1 / d0
    
    return(ans)
}

model.2.cov.func <- function(g, par) {
    
    p.pos <- par[6]
    p0 <- par[5]
    p.neg <- 1 - p0 - p.pos
    
    mu.pos <- par[1] / par[2]
    mu.neg <- - par[3] / par[4]

    a1 <- (mu.neg * p.neg / (p0 + p.neg))^2 * (1 - p.pos * ((1 + g) - 1 / (1 + g)) / g)
    a2 <- 2 * (mu.neg * p.neg / (p0 + p.neg)) * mu.pos * p.pos / (1 + g)
    a3 <- mu.pos * mu.pos * p.pos * g / (1 + g)
    a4 <- p.neg * mu.neg + p.pos * mu.pos

    res <- a1 + a2 + a3 - a4^2
    
    return(res)
}

model.2.est.gamma <- function(cov, par) {
    
    tmp.func <- function(g, cov, par) (cov - model.2.cov.func(g, par))^2
    a <- optimize(f = tmp.func, interval = c(-100, 100), cov = cov, par = par)
    
    return(a)
}

N2G.Spatial.Mixture <- function(data, par.start = c(4, 2, 4, 2, 0.9, 0.05), ksize, ktype = c("2D", "3D"), mask = NULL) {

  ## non-spatial n2g fit ##
  if(is.null(mask)) mask <- data != 0
  data.m <- data[mask == 1]
  fit <- N2G(data.m, par.start)
  ##cat(fit$par)
  p <- fit$par[6]
  
  pos.act.map <- data > fit$lims[1]

  ## lr ##
  lik.ratio <- N2G.Likelihood.Ratio(data, fit$par)
  d <- dim(lik.ratio)
  
  ## neighbourhood (k = 8) ##
  nmat <- matrix(unlist(expand.grid(-1:1, -1:1, 0)), 9, 3)
  nmat <- nmat[-5, ]

  ## number of neighbours ##
  k <- dim(nmat)[1]
  ## possible values of s = number of active voxels in neighbourhood, including central voxel ##
  x <- 0:(k + 1)
  ## number of voxels with s active voxels in the neighbourhood ##
  y <- vector(length = k + 2)

  ## calculate y ##
  for(i in 2:(d[1] - 1))
  {
    for(j in 2:(d[2] - 1))
    {
      for(l in 1:d[3])
      {
        s <- sum(pos.act.map[t(c(i, j, l) + t(nmat))])
        s <- s + pos.act.map[i, j, l]
        y[s + 1] <- y[s + 1] + 1
      }
    }
  }

  ## turn y into a proportion ##
  y <- y / sum(y)

  model.2.func <- function(x, p, gamma, k) {
    ans <- (x == 0) * (1 - (p * ((1 + gamma)^(k + 1) - 1) / (gamma * (1 + gamma)^k)))
    ans <- ans + ((x > 0) * p * gamma^(x - 1) / (1 + gamma)^k)
    
    return(ans)
  }
  
  loglik.func <- function(gamma, x, y, p, k)  {
    -sum(y * log(model.2.func(x, p, gamma, k)))
  }
  
  o <- optim(par = c(gamma = 0.5), fn = loglik.func, method = "L-BFGS-B", lower = c(0), upper = c(Inf), control = list(maxit = 10000, trace = 0), x = x, y = y, k = k, p = p)

  ## parameter estimates ##
  gamma <- o$par[[1]]
  ##cat(gamma)
  
  kt <- switch(ktype[1], "2D" = 2, "3D" = 3)

  ## spatial mixture model ##
  a <- .C("spatial_mixture",
          as.double(aperm(lik.ratio, c(3, 2, 1))),
          as.integer(d),
          as.integer(ksize),
          as.integer(aperm(mask, c(3, 2, 1))),
          as.integer(kt),
          as.double(gamma),
          as.double(p),
          ans = double(prod(d)),
          PACKAGE = "AnalyzeFMRI")
    
  a1 <- array(a$ans, dim = d[3:1])
  a1 <- aperm(a1, c(3, 2, 1))
    
  return(list(p.map = a1, par = fit$par, lims = fit$lims, gamma = gamma, p = p))
}


cov.est <- function(mat, mask, nmat) {

    ## estimate covariance between neighbouring voxels
    
    a <- .C("covariance_est",
            as.double(aperm(mat, c(3, 2, 1))),
            as.integer(dim(mat)),
            as.integer(aperm(mask, c(3, 2, 1))),
            as.integer(t(nmat)),
            as.integer(dim(nmat)),
            ans = double(1),
            PACKAGE = "AnalyzeFMRI")

    return(a$ans)
}




cluster.threshold <- function(x, nmat = NULL, level.thr = 0.5, size.thr) {

  ## thresholds an array at level.thr
  ## calculates the number of contiguous clusters and their sizes
  ## answer is an array in which all voxels that are contained clusters of size greater
  ## than or equal to size.thr are 1, otherwise 0.
  ## nmat is a (Kx3) matrix specifying the neighbourhood system
  ## i.e if a row of nmat is (0, 1, -1) then x[10, 10, 10] and x[10, 11, 9] are neighbours

  if(is.null(nmat)) { ## default is 6 adjacent neighbours
    nmat <- expand.grid(-1:1, -1:1, -1:1)
    nmat <- nmat[c(5, 11, 13, 15, 17, 23), ]
  }   
  
  res <- .C("cluster_mass",
            mat = as.single(aperm(x, c(3, 2, 1))),
            as.integer(dim(x)),
            as.integer(t(nmat)),
            as.integer(dim(nmat)),
            as.single(level.thr),
            num.c = integer(1),
            res.c = single(1000 * 6),
            PACKAGE = "AnalyzeFMRI")

  res.c <- matrix(res$res.c, 1000, 6, byrow = TRUE)[1:res$num.c, ]

  mat1 <- array(res$mat, dim = dim(x)[3:1])
  mat1  <-  aperm(mat1, c(3, 2, 1))

  m <- (res.c[, 5] < size.thr) * (1:res$num.c)
  m <- m[m != 0]

  for(i in 1:length(m)) 
    mat1[mat1 == m[i]] <- 0
  
  mat1 <- 1 * (mat1 > 0)
  
  return(mat1)
}
