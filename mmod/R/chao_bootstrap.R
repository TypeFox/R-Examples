#' Produce bootstrap samples from each subpopulation of a genind object
#'
#' This function produces bootstrap samples from a genind object, with each
#' subpopulation resampled according to its size. Because there are many 
#' statistics that you may wish to calculte from these samples, this function
#' returns a list of genind objects representing bootsrap samples that can then
#' be futher processed (see examples).
#' 
#' You should note, this is a standard (frequentist) approach to quantifying
#' uncertainty - effectively asking "if the population was exactly like our
#' sample, and we repeatedly took samples like this from it, how much would 
#' those samples vary?" The confidence intervals don't include uncertainty 
#' produced from any biases in the way you collected your data. 
#' Additionally, this boostrapping procedure displays a slight upward bias for 
#' some datasets. If you plan or reporting a confidence interval for your 
#' statistic, it is  probably a good idea to subtract the difference between 
#' the point estimate of the statistic and the mean of the boostrap distribution 
#' from the extremes of the interval (as demonstrated in the expample below)
#'
#' @param x genind object (from package adegenet)
#' @param nreps numeric number of bootstrap replicates to perform (default 1000)
#' @export
#' @references Chao, A. et al. (2008). A Two-Stage probabilistic approach to Multiple-Community similarity indices. Biometrics, 64:1178-1186
#' @return A list of genind objects
#' @examples
#'\dontrun{  
#' data(nancycats)
#' obs.D <- D_Jost(nancycats)
#' bs <- chao_bootstrap(nancycats)
#' bs_D <- summarise_bootstrap(bs, D_Jost)
#' bias <- bs.D$summary.global.het[1] - obs.D$global.het
#' bs.D$summary.global.het - bias
#'}
#' @family resample
#' 
#'

chao_bootstrap <- function(x, nreps=1000){
  loc <- x@loc.fac
  pop_sizes <- table(pop(x))
  pop_n <- length(pop_sizes)
  #population allele frequencies
  m <- apply(x@tab,2, function(y) tapply(y, pop(x), mean, na.rm=TRUE) )
  per_pop <- function(i){
    t(do.call(rbind, tapply(m[i,], loc, function(p) rgenotypes(pop_sizes[i],2,p)) ))
  }
  
  one_rep <- function(){
    temp <- x
    temp@tab <- do.call(rbind, sapply(1:pop_n , per_pop, simplify=FALSE))
#    temp@tab <- as.integer(temp@tab)
    return(temp)
  }
  
  return(replicate(nreps, one_rep()))
}
