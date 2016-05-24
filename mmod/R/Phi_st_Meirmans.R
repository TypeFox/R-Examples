#' Calculate Phi_st from a genind object
#'
#' This function calculates Meirmans' corrected version of Phi_st, an Fst 
#' analog produced using the AMOVA framework. Note, the global estimate produced
#' by this function is calculated as the mean distance between individuals
#' across all loci, and this exlcuded individuals with one or more missing 
#' value. 
#'
#' @param x genind object (from package adegenet)
#'
#' @return per.locus Phi_st estimate for each locus
#' @return global Phi_st estimate across all loci 
#' @export
#' @references
#'  Meirmans, PW. (2005), Using the AMOVA framework to estimate a standardized genetic differentiation measure. Evolution 60: 2399-402.
#' @references
#'  Excoffier, L., Smouse, P., Quattro, J. (1992), Analysis of molecular variance inferred from metric distances among DNA haplotypes: application to human mitochondrial DNA restriction data. Genetics 131: 479-91
#' @family diffstat
#' @examples
#' data(nancycats)
#' Phi_st_Meirmans(nancycats[1:26,])



Phi_st_Meirmans <- function(x){


    max_between_dist <- function(pop_name, pop_freqs){
	inter_dists <- sapply(pop_freqs[names(pop_freqs) != pop_name], function(x) x * pop_freqs[pop_name])
	return(sum(inter_dists))
    }
    
	amova_st <- function(dm, dropped){
	
	if(is.null(dropped)){
		pops <- x@pop
	}
	else {
		pops <- droplevels(x@pop[-1 * dropped])
	}

	#values used multiple times
	sq_distance <- dm^2
	pop.freqs <- table(pops)
	n <- dim(dm)[1]
	df <- length(unique(pops)) - 1 
	within_dists <- sapply(unique(pops), function(p) sum(sq_distance[outer(pops==p, pops==p) == 1] ))
	ncoef <- (n - sum(pop.freqs^2)/n)/df

	
	#normal old AMOVA
	SSD_total <- sum(sq_distance/(2 * n))
	SSD_within <- sum(within_dists / (as.numeric(pop.freqs) *2))
	SSD_among <- SSD_total - SSD_within
  
	MSD_total <-SSD_total / (n-1)
	MSD_among <- SSD_among / df
	MSD_within <- SSD_within / ((n -1 ) - df)
	sigma2a <- (MSD_among - MSD_within)/ncoef
	phi <- sigma2a /(sigma2a + MSD_within)
	
	#Now AMOVA to the max
	max_between <- sapply(unique(pops), max_between_dist, pop.freqs)
	max_SSD_total <- sum(within_dists + max_between) / (2*n)
	max_MSD_among <- (max_SSD_total - SSD_within) / df
	sigma2a_prime <- (max_MSD_among - MSD_within)/ncoef
	phi_max = sigma2a_prime /(sigma2a_prime + MSD_within)
	phi_prime <- phi/phi_max
	
	return(phi_prime)
	}

	global <- with(dist.codom(x), amova_st(distances, dropped))
	loci <- sapply(dist.codom(x, global=FALSE), 
	         function(d) with(d, amova_st(distances, dropped)))
	return(list(per.locus=loci, global=global))
}

