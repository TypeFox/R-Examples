#' @name laja
#' @title Macroinvertebrate samples from the Rio Laja of Mexico
#' @description This data set includes macroinvertebrate samples from
#' the Rio Laja, a phylogenetic tree of the taxa and traits that
#' include mean body length and fish feeding preference as in Helmus
#' \emph{et al.} 2013.
#' @docType data
#' @aliases invert.traits invert.tree river.env river.sites
#' @keywords datasets
#' @usage laja
#' @format \code{laja} contains a \code{\link{phylo}} object, a
#' dataframe of sites-by-taxa, a dataframe of sites-by-environment,
#' and a dataframe of traits
#' @references Helmus M., Mercado-Silva N. & Vander Zanden
#' M.J. (2013). Subsidies to predators, apparent competition and the
#' phylogenetic structure of prey communities. Oecologia, 173,
#' 997-1007.
#' @author M.R. Helmus
#' @examples
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits, river.env)
NULL
