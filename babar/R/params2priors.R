UniformPrior = structure(function
### Generate samples from uniform distribution
(u,
### Input scalar/vector of values between 0 and 1
 lowerBound,
### lower bound
 upperBound
### upper bound
 ) {
  return(lowerBound + u*(upperBound - lowerBound))
### Samples from uniform distribution within given bounds
},ex=function(){
    samples <- UniformPrior(runif(1000), -5, 5)
    summary(samples)
})

GaussianPrior = structure(function
### Generate samples from normal distribution
(u,
### Input scalar/vector of values between 0 and 1
 mean,
### mean
 sd
### standard deviation
 ) {
  return(mean + sd*qnorm(u))
### Samples from normal distribution
},ex=function(){
    samples <- GaussianPrior(runif(100), 5, 1)
    summary(samples)
})

LogNormalPrior = structure(function
### Generate samples from the log normal distribution
(u,
### Input scalar/vector of values between 0 and 1
 meanlog,
### Log of mean
 sdlog
### Log of standard deviation
 ) {
  return(exp(meanlog + sdlog*qnorm(u)))## , lower=FALSE)))
### Samples from log normal distribution
},ex=function(){
  samples <- LogNormalPrior(runif(100), 2, 2)
  summary(samples)
})

JeffreysPrior = structure(function
### Generate samples from the Jeffreys prior
(u,
### Input scalar/vector of values between 0 and 1
 log10lowerBound,
### Base 10 logarithm of the lower bound
 log10upperBound
### Base 10 logarithm of the upper bound
 ) { ## So if you have 0.01, pass -2, or if max is 100000000 pass 8
  return(10**(log10lowerBound + u*(log10upperBound - log10lowerBound)))
### Samples from Jeffreys prior
}, ex=function() {
  samples <- JeffreysPrior(runif(100), -2, 8)
  summary(samples)
})

CauchyPrior = structure(function
### Generate samples from the Cauchy prior
(u,
### Input scalar/vector of values between 0 and 1
 location,
### The location parameter (the peak of the distribution)
 scale
### The scale factor (full width at half maximum)
 ) {
  return(location + scale*tan(pi*(u - 0.5)))
### Samples from Cauchy prior
}, ex=function() {
  samples <- CauchyPrior(runif(100), 5, 2)
  summary(samples)
})

ExponentialPrior = structure(function
### Generate samples from the exponential prior
(u,
### Input scalar/vector of values between 0 and 1
 rate
### The rate paramter (such that the mean of the distribution is 1 / rate)
 ) { ## mean = 1/rate
  return(-log(u)/rate)
### Samples from the Exponential prior
}, ex=function() {
  samples <- ExponentialPrior(runif(100), 5)
  summary(samples)
})

WeibullPrior = structure(function
### Generate samples from the Weibull prior
(u,
### Input scalar/vector of values between 0 and 1
 shape,
### The shape parameter of the distriubtion
 scale
### The scale parameter of the distribution
 ) {
  return(scale*(-log(u))**(1/shape))## http://www.johndcook.com/cpp_TR1_random.html#weibull
### Samples from the Weibull prior
}, ex=function() {
  samples <- WeibullPrior(runif(100), 5, 5)
  summary(samples)
})
