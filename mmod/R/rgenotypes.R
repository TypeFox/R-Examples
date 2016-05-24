#' Randomly create genotypes
#'
#' Use the multinomial distribution to randomly create genotpes for individuals
#' for given allele frequences. By default this function returns a matrix of 
#' with alleles in rows and individuals in columns. There is an option to return
#' a genind object representing the same data (see examples).
#' 
#' Used in \code{\link{chao_bootstrap}}, also exported as it may come in handy
#' for other simulations.
#' 
#' @param n integer number of indviduals.
#' @param ploidy integer number of alleles to asign to each individual.
#' @param probs vector of probabilies corresponding to allele frequences.
#' @param genind boolean if TRUE return a genind object
#' @param pop_name charcter Name for population defined in genind object
#' (not required if genind is not TRUE)
#' @param loc_name character name to five locus in genind object
#' @importFrom stats rmultinom
#' @return Either a matrix with individuals in columns, alleles in rows or, if
#' genind is TRUE a genind object for one population and locus.
#' @seealso \code{\link{rmultinom}} which this function wraps.
#' @export
#' @examples
#' 
#' data(nancycats)
#' obs_allele_freqs <- apply(nancycats$tab[,1:16], 2,mean, na.rm=TRUE)
#' rgenotypes(10, 2, obs_allele_freqs)

rgenotypes <- function(n, ploidy, probs, genind=FALSE, pop_name="A", loc_name = "L1"){
 if(all(is.na(probs))){ 
  res <- matrix(NA, ncol=n, nrow=length(probs))
  }
 else
 res <- rmultinom(n, ploidy, probs)
 if(genind){
    res <- t(res)
    colnames(res) <- paste(loc_name, 1:length(probs), sep=".")
    res <- genind(res, rep(pop_name, n, ploidy=ploidy)) 
 }
 return(res)
}






