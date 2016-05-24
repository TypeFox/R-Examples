#' Degree of self belonging of species
#' 
#' Define the degree of self belonging of species.
#' 
#' For the calculation of self-belonging of a set of species the dissimilarities between
#' the species are transformed into similarities and used to define degrees of belonging 
#' to fuzzy sets (Pillar et al. 2009; Pillar & Duarte 2010). Every species among all 
#' species specifies a fuzzy set in relation to all other species, with a certain degree
#'  of belonging. The self-belonging of a given species i expresses its degree of 
#' belonging to the root node of the phylogenetic/functional tree, conditioned to the 
#' similarities between i and all other internal nodes connecting it to the root.
#' 
#' @encoding UTF-8
#' @aliases self.belonging
#' @param dis Matrix containing distance between species.
#' @param standardize Logical argument (TRUE or FALSE) to specify if dis must be standardize
#' in values into range 0 from 1 (Default standardize = TRUE).
#' @return The self-belonging for each species.
#' @author Vanderlei Julio Debastiani <vanderleidebastiani@@yahoo.com.br>
#' @seealso \code{\link{belonging}}
#' @references Pillar, V.D.; Duarte, L.d.S. (2010). A framework for
#' metacommunity analysis of phylogenetic structure. Ecology Letters, 13,
#' 587:596.
#' 
#' Pillar, V.D., Duarte, L.d.S., Sosinski, E.E. & Joner, F. (2009).
#' Discriminating trait-convergence and trait-divergence assembly patterns in
#' ecological community gradients. Journal of Vegetation Science, 20, 334:348.
#' @keywords PCPS
#' @examples
#' 
#' data(flona)
#' self.belonging(flona$phylo)
#' 
#' @export
self.belonging<-function (dis,standardize=TRUE){
	diag.matrix<-diag(belonging(dis,standardize=standardize))
	return(diag.matrix)
}