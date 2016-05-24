#' Bootstrap trees from woodmouse dataset
#'
#' These trees were created using the neighbour-joining and bootstrapping
#' example from the ape documentation.
#'
#'
#' @name woodmiceTrees
#' @docType data
#' @format A multiPhylo object containing 201 trees, each with 15 tips
#' @references Michaux, J. R., Magnanou, E., Paradis, E., Nieberding, C. and
#' Libois, R. (2003) Mitochondrial phylogeography of the Woodmouse (Apodemus
#' sylvaticus) in the Western Palearctic region. \emph{Molecular Ecology}, 12,
#' 685-697
#' @source A set of 15 sequences of the mitochondrial gene cytochrome b of the
#' woodmouse (Apodemus sylvaticus) which is a subset of the data analysed by
#' Michaux et al. (2003). The full data set is available through GenBank
#' (accession numbers AJ511877 to AJ511987)
#' @keywords datasets
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
NULL


#' BEAST analysis of seasonal influenza (A/H3N2)
#'
#' These trees were created using BEAST on hemagglutinin (HA) segments
#' of seasonal influenza A/H3N2 samples collected in New-York city (US)  between 2000 and 2003. This data comes from the influenza BEAST tutorial distributed at:
#' http://beast.bio.ed.ac.uk/tutorials
#'
#' Only the first 200 trees (out of 10,000) were retained.
#'
#' @name fluTrees
#' @docType data
#' @format A multiPhylo object containing 200 trees, each with 165 tips
#' @references http://beast.bio.ed.ac.uk/tutorials
#' @source http://beast.bio.ed.ac.uk/tutorials
#' @keywords datasets
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
NULL

#' BEAST analysis of Dengue fever
#'
#' These trees were created using one of the \code{xml} files provided with the original BEAST paper by Drummond and Rambaut (2007).
#' They provide an example of 17 dengue virus serotype 4 sequences from Lanciotti et al. (1997) and \code{xml} files with varying priors for model and clock rate.
#' Here we include a random sample of 500 of the trees (from the second half of the posterior) produced using BEAST v1.8 with the standard GTR + Gamma + I substitution model with uncorrelated lognormal-distributed relaxed molecular clock (file 4).
#'
#' @name DengueTrees
#' @docType data
#' @format A multiPhylo object containing 500 trees, each with 17 tips
#' @references Drummond, A. J., and Rambaut, A. (2007) 
#' BEAST: Bayesian evolutionary analysis by sampling trees.
#' \emph{BMC Evolutionary Biology}, 7(1), 214.
#'
#' Lanciotti, R. S., Gubler, D. J., and Trent, D. W. (1997)
#' Molecular evolution and phylogeny of dengue-4 viruses.
#' \emph{Journal of General Virology}, 78(9), 2279-2286.
#' @source http://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-7-214
#' @keywords datasets
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
NULL