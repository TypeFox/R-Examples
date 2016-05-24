.generateTransform <- function(data) {
  mu = mean(data)
  sigma = sd(data)

  transform <- function(params) {
    tParams = numeric(length=length(params))
    tParams[1] = GaussianPrior(params[1], mu, sigma)
    tParams[2] = UniformPrior(params[2], 0, 2 * sigma)

    return(tParams)
  }

  return(transform)
  
}

.llfMaker <- function(data, transformParams) {
  llf <- function(params) {

    tParams = transformParams(params)
    
    mean = tParams[1]
    sigma = tParams[2]

    n <- length(data)

    ll <- -(n / 2) * log(2 * pi) - (n / 2) * log(sigma ** 2) - (1 / (2 * sigma ** 2)) * sum((data - mean) ** 2)

    return(ll)
  }

  return(llf)
}

compareDistributions <- structure(function(
### Compare two sets of normally distributed samples using nested sampling,
### to determine whether they have the same mean and variance.
  data.first,
### Samples from the first distriubtion, as a vector of normally distributed values
  data.second
### Samples from the second distribution as a vector of normally distributed values
  ) {

  transformParams.first <- .generateTransform(data.first)
  transformParams.second <- .generateTransform(data.second)
  
  llf.first <- .llfMaker(data.first, transformParams.first)
  llf.second <- .llfMaker(data.second, transformParams.second)

  data.combined <- c(data.first, data.second)
  transformParams.combined <- .generateTransform(data.combined)
  
  llf.combined <- .llfMaker(data.combined, transformParams.combined)

  prior.size <- 25
  tol <- 0.5
  
  ns.first <- nestedSampling(llf.first, 2, prior.size, transformParams.first, tolerance=tol)
  ns.second <- nestedSampling(llf.second, 2, prior.size, transformParams.second, tolerance=tol)
  ns.combined <- nestedSampling(llf.combined, 2, prior.size, transformParams.combined, tolerance=tol)

  evidence.different <- ns.first$logevidence + ns.second$logevidence
  evidence.same <- ns.combined$logevidence

  return(evidence.same - evidence.different)
  ### Bayes factor for the hypothesis that the distributions have the same mean and variance
  ### versus the hypothesis that they have different means and variances
},ex=function() {
    data.a <- rnorm(10, 1, 1)
    data.b <- rnorm(10, 5, 1)
    compareDistributions(data.a, data.b)
})
