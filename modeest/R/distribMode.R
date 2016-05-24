# Author : P. Poncet

########################################################################
# Mode (most likely value) of some continuous and discrete distributions
########################################################################

.distribList <- c("beta",
                  "cauchy",
                  "chisq",
                  "exp",
                  #"extreme",
                  "f",
                  "fisk", 
                  "frechet",
                  "gamma",
                  "norm",
                  "gh",
                  "gev",
                  "gompertz", 
                  "gpd",
                  "gumbel",
                  "hyp",
                  "kumar",
                  "laplace",
                  "koenker", 
                  "logis",
                  #"log_logistic",
                  "lnorm",
                  #"maxwell_boltzmann", 
                  "nig",
                  "paralogistic", 
                  "pareto",
                  "rayleigh",
                  "rweibull",
                  "stable",
                  #"symstb",
                  "t",
                  "unif",
                  "weibull",
                  "bern",
                  "binom",
                  "geom",
                  "hyper",
                  #"logarithmic",
                  #"multinomial",
                  "nbinom",
                  "pois")
                  #"signrank",
                  #"wilcoxon")

         
#######################################################
## Continuous distributions
#######################################################

#------------------------------------------------------

### Beta distribution

betaMode <-
function(shape1,
         shape2,
         ncp = 0)
{
  if (ncp == 0) {
    M <- (shape1-1)/(shape1+shape2-2)
  } else {
    warning("still to be done. 'NA' is returned")
    M <- NA
  }
  return(M)
}

.mlv.beta <-
function(...)
{
  M <- betaMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = 1 - 2*pbeta(M, ...),
                                  x = "beta",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))  
}


#------------------------------------------------------

### Cauchy distribution

cauchyMode <-
function(location = 0,
         ...)
{
  return(location)
}

.mlv.cauchy <-
function(...)
{
  M <- cauchyMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = 1 - 2*pcauchy(M, ...),
                                  x = "cauchy",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Chernoff distribution

#chernMode <- 
#function(location = 0,
#         ...)
#{
#  return(location)
#} 

#.mlv.chernoff   ...


#------------------------------------------------------


### Chi-square distribution

chisqMode <-
function(df,
         ncp=0)
{
  if (ncp == 0 & df == 0) {
    M <- 0
  } else if (ncp == 0 & df %in% c(1,2)) {
    warning("the density is not continuous at the mode.")
    M <- 0
  } else if (ncp == 0 & df >= 3) {
    M <- df - 2
  } else {
    warning("still to be done. 'NA' is returned")
    M <- NA
  }
  return(M)
}

.mlv.chisq <-
function(...)
{
  M <- chisqMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = 1 - 2*pchisq(M, ...),
                                  x = "chisquare",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Exponential distribution
expMode <-
function(...)
{
  return(0)
}

.mlv.exp <-
function(...)
{
  M <- expMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = 1 - 2*pexp(M, ...),
                                  x = "exponential",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Extreme distribution (related to package 'survival')

#extremeMode <-
#function(loc=0,
#         scale=1,
#         shape=0)
#{
#  NA
#}

#.mlv.extreme <-
#function(...)
#{
#  M <- extremeMode(...)
#  return(invisible(structure(list(M = M,
#                                  skewness = 1 - 2*pextreme(M, ...),
#                                  x = "extreme",
#                                  method = "continuous",
#                                  call = match.call()),
#                             class = "mlv")))
#}


#------------------------------------------------------

### F distribution

fMode <-
function(df1,
         df2)
{ 
  if (df1 > 2) {
    M <- (1-2/df1)*(df2/(2+df2))
  } else {
    warning("still to be done. 'NA' is returned")
    M <- NA
  }
  return(M)
}

.mlv.f <-
function(...)
{
  M <- fMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = 1 - 2*pf(M, ...),
                                  x = "F",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Fisk distribution

fiskMode <- 
function(shape1.a, 
         scale = 1)
{
  return(scale*((shape1.a-1)/(shape1.a+1))^(1/shape1.a))
}

.mlv.fisk <-
function(...)
{
  #require(VGAM)
  if (!"package:VGAM" %in% search()) {
    warning("package 'VGAM' should be loaded.")
  }
  M <- fiskMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = try(1 - 2*pfisk(M, ...), silent = TRUE),
                                  x = "fisk",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Fréchet distribution

frechetMode <-
function(loc=0,
         scale=1,
         shape=1,
         ...)
{
  #return(loc + (scale/abs(shape))*((1+abs(shape))^(-abs(shape))-1))
  return(loc + scale*(shape/(1+shape))^(1/shape))
}

.mlv.frechet <-
function(...)
{
  #require(evd)
  if (!"package:evd" %in% search()) {
    warning("package 'evd' should be loaded.")
  }
  M <- frechetMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = try(1 - 2*pfrechet(M, ...), silent = TRUE),
                                  x = "frechet",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Gamma distribution

gammaMode <-
function(shape,
         rate = 1,
         scale = 1/rate)
{
  return(scale*(shape-1))
}

.mlv.gamma <-
function(...)
{
  M <- gammaMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = 1 - 2*pgamma(M, ...),
                                  x = "gamma",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Gaussian (normal) distribution

normMode <-
function(mean = 0,
         ...)
{
  return(mean)
}

.mlv.norm <-
function(...)
{
  M <- normMode(...) 
  return(invisible(structure(list(M = M,
                                  skewness = 1 - 2*pnorm(M, ...),
                                  x = "gaussian",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}



#------------------------------------------------------

### Generalized extreme value distribution

#! Changement
#gevMode <-
#function(loc=1,
#         scale=1,
#         shape=1,
#         ...)
#{
#  if (shape==0) {
#    M <- loc
#  } else {
#    M <- loc + (scale/shape)*(max(0,(1+shape))^(-shape)-1)
#  }
#  return(M)  
#}

gevMode <-
function(loc=0,
         scale=1,
         shape=0,
         ...)
{
  k <- pmax(0,(1+shape))^(-shape)-1 #! il y avait une erreur ici : je mets 'pmax' au lieu de 'max'
  shape[shape==0] <- Inf
  M <- loc + (scale/shape)*k
  return(M)  
}


.mlv.gev <-
function(...)
{
  #require(evd)
  if (!"package:evd" %in% search()) {
    warning("package 'evd' should be loaded.")
  }
  M <- gevMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = try(1 - 2*pgev(M, ...), silent = TRUE),
                                  x = "gev",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))  
}


#------------------------------------------------------

### Generalized hyperbolic distribution

## The following function is taken from package 'fBasics'
ghMode <-
function(alpha = 1,
         beta = 0,
         delta = 1,
         mu = 0,
         lambda = -1/2,
         ...)
{
  min <- qgh(0.01, alpha, beta, delta, mu, lambda)
  max <- qgh(0.99, alpha, beta, delta, mu, lambda)
  M <- optimize(f = dgh, interval = c(min, max), alpha = alpha, 
      beta = beta, delta = delta, mu = mu, lambda = lambda, 
      maximum = TRUE, tol = .Machine$double.eps)$maximum
  return(M)
}

.mlv.gh <-
function(...)
{
  #require(fBasics)
  if (!"package:fBasics" %in% search()) {
    warning("package 'fBasics' should be loaded.")
  }
  M <- ghMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = try(1 - 2*pgh(M, ...), silent = TRUE),
                                  x = "generalized_hyperbolic",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))  
}


#------------------------------------------------------

### Generalized Pareto distribution

#! Attention au cas weibull ?
gpdMode <-
function(loc=0,
         scale=1,
         shape=0,
         ...)
{
  if (shape==-1) {
    warning("all values between 'loc' and 'loc+scale' are modes. Only the mean value is returned")
    M <- loc + scale/2  
  } else if (-2-1/shape > 0) {
    M <- loc - scale/shape
  } else {
    warning("the density is not continuous at the mode.")
    M <- loc
  }
  return(M)
}

.mlv.gpd <-
function(...)
{
  #require(evd)
  if (!"package:evd" %in% search()) {
    warning("package 'evd' should be loaded.")
  }
  M <- gpdMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = try(1 - 2*pgpd(M, ...), silent = TRUE),
                                  x = "gpd",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Gompertz distribution

gompertzMode <-
function(shape,
         scale=1)
{
  if (shape < scale) {
    M <- log(scale/shape)/scale
  } else {
    M <- 0
  }
  return(M)
}

.mlv.gompertz <-
function(...)
{
  #require(VGAM)
  if (!"package:VGAM" %in% search()) {
    warning("package 'VGAM' should be loaded.")
  }
  M <- gompertzMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = try(1 - 2*pgompertz(M, ...), silent = TRUE),
                                  x = "gompertz",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}

#------------------------------------------------------

### Gumbel distribution

gumbelMode <-
function(loc=0,
         ...)
{
  return(loc)  
}

.mlv.gumbel <-
function(...)
{
  #require(evd)
  if (!"package:evd" %in% search()) {
    warning("package 'evd' should be loaded.")
  }
  M <- gumbelMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = try(1 - 2*pgumbel(M, ...), silent = TRUE),
                                  x = "gumbel",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))  
}


#------------------------------------------------------

### Hyperbolic distribution

## The 5 following functions are taken from package 'fBasics'

.hyp1Mode <-
function(alpha = 1,
         beta = 0,
         delta = 1,
         mu = 0) 
{
  return(mu + delta * beta/sqrt(alpha^2 - beta^2))
}

.hyp2Mode <-
function(zeta = 1,
         rho = 0,
         delta = 1,
         mu = 0) 
{
  alpha <- zeta / ( delta * sqrt(1 - rho*rho) )
  return(hypMode(alpha, alpha * rho, delta, mu))
}

.hyp3Mode <-
function(xi = 1/sqrt(2),
         chi = 0,
         delta = 1,
         mu = 0) 
{
  rho <- chi/xi
  zeta <- 1/xi^2 - 1
  alpha <- zeta/(delta * sqrt(1 - rho * rho))
  beta <- alpha * rho
  return(hypMode(alpha, beta, delta, mu))
}

.hyp4Mode <-
function(a.bar = 1,
         b.bar = 0,
         delta = 1,
         mu = 0) 
{
  return(hypMode(a.bar/delta, b.bar/delta, delta, mu))
}

hypMode <-
function(alpha = 1,
         beta = 0,
         delta = 1,
         mu = 0,
         pm = c(1, 2, 3, 4)) 
{
  return(eval(call(paste(".hyp", pm[1], "Mode", sep = ""), alpha, beta, delta, mu)))
}

.mlv.hyp <-
function(...)
{
  #require(fBasics)
  if (!"package:fBasics" %in% search()) {
    warning("package 'fBasics' should be loaded.")
  }
  M <- hypMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = try(1 - 2*phyp(M, ...), silent = TRUE),
                                  x = "hyperbolic",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))  
}


#------------------------------------------------------

### Koenker distribution 

koenkerMode <- 
function(location = 0, 
         ...)
{
  return(location)
}

.mlv.koenker <- 
function(...)
{
  #require(VGAM)
  if (!"package:VGAM" %in% search()) {
    warning("package 'VGAM' should be loaded.")
  }
  M <- koenkerMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = try(1 - 2*pkoenker(M, ...), silent = TRUE),
                                  x = "koenker",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Kumaraswamy distribution

kumarMode <-
function(shape1,
         shape2)
{
  return((shape1-1)/(shape1*shape2 - 1)^(1/shape1))
}

.mlv.kumar <-
function(...)
{
  #require(VGAM)
  if (!"package:VGAM" %in% search()) {
    warning("package 'VGAM' should be loaded.")
  }
  M <- kumarMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = try(1 - 2*pkumar(M, ...), silent = TRUE),
                                  x = "kumaraswamy",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Laplace distribution

laplaceMode <-
function(location = 0,
         ...)
{
  return(location)
}

.mlv.laplace <-
function(...)
{
  #require(VGAM)
  if (!"package:VGAM" %in% search()) {
    warning("package 'VGAM' should be loaded.")
  }
  M <- laplaceMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = try(1 - 2*plaplace(M, ...), silent = TRUE),
                                  x = "laplace",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Logistic distribution

logisMode <-
function(location = 0,
         ...)
{
  return(location)
}

.mlv.logis <-
function(...)
{
  M <- logisMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = 1 - 2*plogis(M, ...),
                                  x = "logistic",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Log-logistic distribution (related to package 'survival')

#loglogisMode <-
#function()
#{
#  NA
#}

#.mlv.log_logistic
#function(...)
#{
#  M <- loglogisMode(...)
#  return(invisible(structure(list(M = M,
#                                  skewness = 1 - 2*ploglogis(M, ...),
#                                  x = "log_logistic",
#                                  method = "continuous",
#                                  call = match.call()),
#                             class = "mlv")))
#}


#------------------------------------------------------

### Lognormal distribution

lnormMode <-
function(meanlog = 0,
         sdlog = 1)
{
  return(exp(meanlog - sdlog^2))
}

.mlv.lnorm <-
function(...)
{
  M <- lnormMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = 1 - 2*plnorm(M, ...),
                                  x = "lognormal",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}



#------------------------------------------------------

### Maxwell-Boltzmann distribution

#mbMode <-
#function(scale,
#         ...)
#{
#  return(scale*sqrt(2))
#}

#.mlv.maxwell_boltzmann <-
#function(...)
#{
#  M <-  mbMode(...)
#  return(invisible(structure(list(M = M,
#                                  skewness = 1 - 2*NA,
#                                  x = "maxwell_boltzmann",
#                                  method = "continuous",
#                                  call = match.call()),
#                             class = "mlv")))
#}


#------------------------------------------------------

### Normal Inverse Gaussian distribution

## The following function is taken from package 'fBasics'

nigMode <-
function(alpha = 1,
         beta = 0,
         delta = 1,
         mu = 0,
         ...)
{
  min <- qnig(0.01, alpha, beta, delta, mu)
  max <- qnig(0.99, alpha, beta, delta, mu)
  M <- optimize(f = dnig, interval = c(min, max), alpha = alpha, 
      beta = beta, delta = delta, mu = mu, maximum = TRUE, 
      tol = .Machine$double.eps)$maximum
  return(M)
}

.mlv.nig <-
function(...)
{
  #require(fBasics)
  if (!"package:fBasics" %in% search()) {
    warning("package 'fBasics' should be loaded.")
  }
  M <- nigMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = try(1 - 2*pnig(M, ...), silent = TRUE),
                                  x = "normal_inverse",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Paralogistic distribution

paralogisticMode <- 
function(shape1.a, 
         scale = 1)
{
  if (shape1.a <= 1) {
    warning("the density is not continuous at the mode.")
    return(0)
  } else {
    return(scale*((shape1.a-1)/(shape1.a^2+1))^(1/shape1.a))
  }
}

.mlv.paralogistic <-
function(...)
{
  #require(VGAM)
  if (!"package:VGAM" %in% search()) {
    warning("package 'VGAM' should be loaded.")
  }
  M <- paralogisticMode(...) 
  return(invisible(structure(list(M = M,
                                  skewness = try(1 - 2*pparalogistic(M, ...), silent = TRUE),
                                  x = "paralogistic",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Pareto distribution

paretoMode <-
function(location,
         ...)
{
  warning("the density is not continuous at the mode.")  
  return(location)
}

.mlv.pareto <-
function(...)
{
  #require(VGAM)
  if (!"package:VGAM" %in% search()) {
    warning("package 'VGAM' should be loaded.")
  }
  M <- paretoMode(...) 
  return(invisible(structure(list(M = M,
                                  skewness = try(1 - 2*ppareto(M, ...), silent = TRUE),
                                  x = "pareto",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Rayleigh distribution       

rayleighMode <-
function(scale = 1)
{
  return(scale)
}

.mlv.rayleigh <-
function(...)
{
  #require(VGAM)
  if (!"package:VGAM" %in% search()) {
    warning("package 'VGAM' should be loaded.")
  }
  M <- rayleighMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = try(1 - 2*prayleigh(M, ...), silent = TRUE),
                                  x = "rayleigh",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Stable distribution

## The following function is taken from package 'stabledist'

stableMode <-
function(alpha,
         beta,
         gamma = 1,
         delta = 0,
         pm = 0, 
         ...) 
{
  beta.max <- 1 - 1e-11
  tol <- .Machine$double.eps^0.25
  if (gamma == 1 & delta == 0 & pm == 0) {
    stopifnot(0 < alpha, alpha <= 2, length(alpha) == 1, -1 <= beta, beta <= 1, length(beta) == 1, length(beta.max) == 1)
    if (alpha * beta == 0) {
      M <- 0
    }
    if (beta > beta.max) {
      beta <- beta.max
    }
    M <- optimize(dstable, interval = c(-0.7, 0) * sign(beta), alpha = alpha, 
        beta = beta, pm = 0, maximum = TRUE, tol = tol)$maximum    
    return(M)
  } else {
    warning("still to be done. 'NA' is returned")
    return(NA)
  }
}

.mlv.stable <-
function(...)
{
  #require(stabledist)
  if (!"package:stabledist" %in% search()) {
    warning("package 'stabledist' should be loaded.")
  }
  M <- stableMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = try(1 - 2*pstable(M, ...), silent = TRUE),
                                  x = "stable",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Symmetric stable distribution

#symstbMode <-
#function(...)
#{
#  return(0)
#}
#
#.mlv.symstb <-
#function(...)
#{
#  M <- symstbMode(...)
#  return(invisible(structure(list(M = M,
#                                  skewness = 1 - 2*psymstb(M, ...),
#                                  x = "symmetric_stable",
#                                  method = "continuous",
#                                  call = match.call()),
#                             class = "mlv")))
#}


#------------------------------------------------------

### Weibull distribution (in the context of extreme value theory)

rweibullMode <-
function(loc = 0,
         scale = 1,
         shape = 1,
         ...)
{
  #return(loc + (scale/-abs(shape))*((1-abs(shape))^(abs(shape))-1))
  if (shape < 1) {
    M <- loc
  } else {
    M <- loc - scale*((shape-1)/shape)^(1/shape)
  }
}

.mlv.rweibull <-
function(...)
{
  #require(evd)
  if (!"package:evd" %in% search()) {
    warning("package 'evd' should be loaded.")
  }
  M <- rweibullMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = try(1 - 2*prweibull(M, ...), silent = TRUE),
                                  x = "rweibull",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Student distribution

tMode <-
function(df,
         ncp = 0)
{
  if (ncp == 0) {
    M <- 0
  } else {
    warning("still to be done. 'NA' is returned")
    M <- NA
  }
  return(M)
}

.mlv.t <-
function(...)
{
  M <- tMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = 1 - 2*pt(M, ...),
                                  x = "student",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}



#------------------------------------------------------

### Uniform distribution

unifMode <-
function(min = 0,
         max = 1)
{
  warning("all values between 'min' and 'max' are modes. Only the mean value is returned")
  return((min+max)/2)
}

.mlv.unif <-
function(...)
{
  M <- unifMode(...)
  return(invisible(structure(list(M = M,
                                skewness = 1 - 2*punif(M, ...),
                                x = "uniform",
                                method = "continuous",
                                call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Weibull distribution

weibullMode <-
function(shape,
         scale = 1,
         ...)
{
  return(scale*(1-1/shape)^(1/shape))
}

.mlv.weibull <-
function(...)
{
  M <- weibullMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = 1 - 2*pweibull(M, ...),
                                  x = "weibull",
                                  method = "continuous",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------




#######################################################
## Discrete distributions
#######################################################

#------------------------------------------------------

### Bernoulli distribution

bernMode <-
function(prob)
{
  if (prob > 1 || prob < 0) return(NaN)
  q <- 1 - prob
  if (q > prob) {
    return(0)
  } else {
    if (q == prob) {
      return(c(0,1))
    } else {
      if (q < prob) return(1)
    }
  }
}

.mlv.bern <-
function(...)
{
  M <- bernMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = 1 - 2*pbinom(M - 1, size = 1,...) - dbinom(M, size = 1,...),
                                  x = "bernoulli",
                                  method = "discrete",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Binomial distribution

binomMode <-
function(size,
         prob)
{
  
  if (prob > 1 || prob < 0) return(NaN)
  if (prob == 0) {
    return(0)
  } else {
    if (prob == 1) {
      return(size)
    } else {
      x <- ceiling((size+1)*prob - 1)
      if (x == (size+1)*prob - 1) {
        return(c(x, x+1))
      } else {
        return(x)
      }      
    }
  }     
}

.mlv.binom <-
function(...)
{
  M <- binomMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = 1 - 2*pbinom(M - 1, ...) - dbinom(M, ...),
                                  x = "binomial",
                                  method = "discrete",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Geometric distribution

geomMode <-
function(...)
{
  return(1)
}

.mlv.geom <-
function(...)
{
  M <- geomMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = 1 - 2*pgeom(M - 1, ...) - dgeom(M, ...),
                                  x = "geometric",
                                  method = "discrete",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Hypergeometric distribution

# value of the mode seen on http://www.math.uah.edu/stat/urn/Hypergeometric.html

hyperMode <-
function(m,
         n,
         k,
         ...)
{
  lambda <- (m+1)*(k+1)/(m+n+1)  
  if (lambda == 0) {
    return(0)
  } else {
    x <- floor(lambda)
    if (lambda == x) {
      return(c(lambda-1,lambda))
    } else {
      return(x)
    }
  }
}

.mlv.hyper <-
function(...)
{
  M <- hyperMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = 1 - 2*phyper(M - 1, ...) - dhyper(M, ...),
                                  x = "hypergeometric",
                                  method = "discrete",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Logarithmic distribution

#logMode <-
#function(...)
#{
#  return(1)
#}

#.mlv.logarithmic <-
#function(...)
#{
#  M <- logMode(...)
#  return(invisible(structure(list(M = M,
#                                  skewness = 1 - 2*NA - NA,
#                                  x = "logarithmic",
#                                  method = "discrete",
#                                  call = match.call()),
#                             class = "mlv")))
#}


#------------------------------------------------------

### Multinomial distribution

#multinomMode <-
#function(size,
#         prob)
#{
#  if (prob > 1 || prob < 0) return(NaN)
#  warning("Still to be done. 'NA' is returned.")
#  return(NA)
#}

#.mlv.multinomial <-
#function(...)
#{
#  M <- multinomMode(...)
#  return(invisible(structure(list(M = M,
#                                  skewness = 1 - 2*NA - NA,
#                                  x = "multinomial",
#                                  method = "discrete",
#                                  call = match.call()),
#                             class = "mlv")))
#}


#------------------------------------------------------

### Negative binomial distribution

nbinomMode <-
function(size,
         prob,
         mu)
{
  if (prob > 1 || prob < 0) return(NaN)
  if (!missing(mu)) {
    prob <- size/(size+mu)    
  }
  if (size <= 1) {
    return(0)
  } else {
    return(floor((size-1)*(1-prob)/prob))
  }
}

.mlv.nbinom <-
function(...)
{
  M <- nbinomMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = 1 - 2*pnbinom(M - 1, ...) - dnbinom(M, ...),
                                  x = "negative_binomial",
                                  method = "discrete",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Poisson distribution

poisMode <-
function(lambda)
{
  if (lambda < 0) return(NaN)
  if (lambda == 0) return(0)
  else {
    x <- floor(lambda)
    if (lambda == x) {
      return(c(lambda-1,lambda))
    } else {
      return(x)
    }
  }
}

.mlv.pois <-
function(...)
{
  M <- poisMode(...)
  return(invisible(structure(list(M = M,
                                  skewness = 1 - 2*ppois(M - 1, ...) - dpois(M, ...),
                                  x = "poisson",
                                  method = "discrete",
                                  call = match.call()),
                             class = "mlv")))
}


#------------------------------------------------------

### Signrank distribution

#signrankMode <-
#function(n,
#         ...)
#{
#  warning("Still to be done. 'NA' is returned.")
#  return(NA)
#}

#.mlv.signrank <-
#function(...)
#{
#  M <- signrankMode(...)
#  return(invisible(structure(list(M = M,
#                                  skewness = 1 - 2*NA - NA,
#                                  x = "signrank",
#                                  method = "discrete",
#                                  call = match.call()),
#                             class = "mlv")))
#}


#------------------------------------------------------

### Wilcoxon distribution

#wilcoxMode <-
#function(m,
#         n,
#         ...)
#{
#  warning("Still to be done. 'NA' is returned.")
#  return(NA)
#}

#.mlv.wilcox <-
#function(...)
#{
#  M <- wilcoxMode(...)
#  return(invisible(structure(list(M = M,
#                                  skewness = 1 - 2*NA - NA,
#                                  x = "wilcoxon",
#                                  method = "discrete",
#                                  call = match.call()),
#                             class = "mlv")))
#}

