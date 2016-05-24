#' Simulation of distance sampling data via mixture models

#' Allows one to simulate line transect distance sampling data using 
#' a mixture of half-normal detection functions.
#'
#' @param n number of samples to generate
#' @param sigma vector of scale parameters
#' @param mix.prop vector of mixture proportions (same length as sigma)
#' @param width truncation
#' @param means vector of means (used to generate wacky, non-monotonic data)
#' @return distances a vector of distances
#'
#' @author David Lawrence Miller
#'
#' @note At the moment this is TOTALLY UNSUPPORTED!
#'       Please don't use it for anything important!

sim.mix<-function(n, sigma,mix.prop,width,means=0){

   # TODO:
   #  * covariates
   #  * From SNW: "All following could be made more efficient by not 
   #               generating so many more than needed..."


   # code from Simon Wood for direct simulation of mixtures

   ## detection function is of form \sum_i p_i exp(-d^2/(2 sigma_i^2))
   ## truncated at width, say. This is equivalent to 
   ## \sum_i \alpha_i /\sqrt{2 \pi sigma_i^2} exp(-d^2/(2 \sigma_i^2))
   ## i.e. to a mixture of half normals with mixture weights proportional to
   ## \alpha_i = p_i \sqrt{2 \pi sigma_i^2}
   ## the truncated version of the latter is easy to simulate from.

   
   
   ## mixture method
   alpha <- mix.prop*(sqrt(2*pi)*sigma) ## transform to mixture weights
   alpha <- alpha/(sum(alpha)) ## normalize
   ii <- sample(1:length(mix.prop),size=2*n,replace=TRUE,prob=alpha) ## sample mixture index
   d <- abs(rnorm(ii,means,sigma[ii])) ## sample from components
   d <- d[d<=width]   ## reject beyond truncation
   d <- d[1:n] ## ditch the excess

   return(d)
}

