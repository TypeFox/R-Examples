# calculate the bandwidth using 
# options to use:
# - bw.ucv( Y )
# - bw.bcv( Y )
# - bw.SJ( Y ) ## maybe the best
# - bw.nrd0( Y ) ## adjusted Silverman

# rule of thumb
bandwidth_rule <- function( dimX, dimY, Nsamples ) {
  return( Nsamples ^ ( -1 / (3 * (dimX + dimY + 1)) ) )
}

## Sampling from 1D mixture of Gaussians without bound (Badin and Simar 2003, page 14)
## n  - number of samples
## Y  - the means of the Gaussians
## bw - the bandwidth of the Gaussians
sampling_from_1d_gaussian_mixture <- function( n, Y, bw ) {
  # select n random indices
  ind   = sample.int( length(Y), n, replace=TRUE )
  # now select the means according to the indeces
  means = Y[ind]
  # generating samples from normal distributions:
  samples = rnorm(n, mean = means, sd = matrix(bw,n,1) )
  # returning samples:
  return( samples )
}

## Reflection sampling method (Silverman, 1986)
## Samples from the range (-Inf, 1]
## n     - number of samples
## theta - the means of the Gaussians (not above 1)
## bw    - the bandwidth of the Gaussians
## tvar  - variance of the original theta samples
## tmean - mean of the original theta samples
sampling_with_reflection <- function( n, theta, bw, tvar, tmean ) {
  samples = sampling_from_1d_gaussian_mixture( n, theta, bw )
  # correct for the samples that are bigger than 1:
  samples[samples > 1.0] = 2 - samples[samples > 1.0]
  samples = tmean + 1/sqrt(1 + bw^2 / tvar) * ( samples - tmean )
  return(samples)
}

## Reflection sampling method that samples in the range [1, Inf)
sampling_delta_with_reflection <- function( n, delta, bw, dvar, dmean ) {
  samples = sampling_from_1d_gaussian_mixture( n, delta, bw )
  # correct samples outside the [1, Inf):
  samples[samples < 1.0] = 2 - samples[samples < 1.0]
  samples = dmean + 1/sqrt(1 + bw^2 / dvar) * (samples - dmean)
  return(samples)
}

sampling_logtheta_with_reflection <- function( n, ltheta, lbw, ltvar, ltmean ) {
  samples = sampling_from_1d_gaussian_mixture( n, ltheta, lbw )
  # correct for the samples that are bigger than 0:
  samples[samples > 0.0] = - samples[samples > 0.0]
  samples = ltmean + 1/sqrt(1 + lbw^2 / ltvar) * ( samples - ltmean )
  return(samples)
}
