#' Phylogenetics for the Environmental Sciences
#'
#' Analysis and manipulation of eco-phylogenetic datasets containing
#' species phylogeny, species traits, community composition, and
#' environmental data. Provide the \code{\link{comparative.comm}}
#' object to ease data manipulation, and wrappers for common community
#' phylogenetic indices grouped according to Pearse et al. 2014:
#' \code{\link{pez.shape}}, \code{\link{pez.evenness}},
#' \code{\link{pez.dispersion}}, and
#' \code{\link{pez.dissimilarity}}. Implementation of Cavender-Bares
#' et al. (2004) correlation of phylogenetic and ecological matrices
#' (\code{\link{fingerprint.regression}}). Simulation of null
#' assemblages, traits, and phylogenies (\code{\link{scape}},
#' \code{\link{sim.meta.comm}}).
#' 
#' @docType package
#' @name pez
#' @aliases pez package-pez pez-package
#' @references Pearse W.D., Purvis A., Cavender-Bares J. & Helmus M.R. (2014). Metrics and Models of Community Phylogenetics. In: Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology. Springer Berlin Heidelberg, pp. 451-464.
#' @references Cavender-Bares J., Ackerly D.D., Baum D.A. & Bazzaz F.A. (2004) Phylogenetic overdispersion in Floridian oak communities. The Americant Naturalist 163(6): 823--843.
#' @examples
#' require(pez)
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits, river.env)
#' pez.shape(data)
NULL

