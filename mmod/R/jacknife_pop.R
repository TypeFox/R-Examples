#' Create jacknife samples of a genind object by population
#'
#' Makes a series of jacknife samples across populations from a genind object.
#' This function returns a list of genind objects that can then be further 
#' processed (see examples below).
#'
#' @param x genind object (from package adegenet)
#' @param sample_frac fraction of pops to sample in each replication (default 0.5)
#' @param nreps number of jacknife replicates to run (default 1000)
#' @export
#' @return  a list of genind objects to be further processed
#' @examples
#'\dontrun{  
#' data(nancycats)
#' obs <- diff_stats(nancycats)
#' jn <- jacknife_populations(nancycats)
#' jn.D <- summarise_bootstrap(jn, D_Jost)
#' }
#' @family resample

jacknife_populations <- function(x, sample_frac=0.5, nreps=1000){
 on.exit(cat("\n")) 
 rep <- function(i,d){
  if( i %% 50 == 0){
    cat("\r", paste(i, "of", nreps, "reps completed"))
  }
  to.sample <- sample(d$pop.names, length(d$pop.names) * sample_frac)
  return(d[d$pop.names %in% to.sample,])
 }
 return(sapply(1:nreps, rep, x))
}
  
