#' Maximum likelihood tree for 13 \emph{Ligophorus} species
#'
#' This data set contains the maximum likelihood (Felsenstein, 1981) phylogenetic tree of 13 \emph{Ligophorus} species, 
#' inferred using concatenated 18S rRNA, ITS1 rRNA and 28S rRNA DNA sequences (10 000 bootstrap replicates). The best DNA substitution 
#' model was GTR+G. Multiple sequence alignment was done using online MAFFT server (Version 7; Katoh, 2013; Katoh & Standley, 2013) 
#' with default parameters. Maximum likelihood tree construction was done using IQ-TREE (Minh et al., 201; Nguyen et al., 2015). 
#' @docType data
#' @usage data(ligotree)
#' @format A \code{phylo} object, so the \code{ape} package (Paradis et al., 2004) needs to be installed
#' @details The phylogenetic tree comes with bootstrap support for the internal nodes.
#' @keywords datasets
#' @source Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#' 
#' Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Data from: Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. Dryad Digital Repository. http://dx.doi.org/10.5061/dryad.50sg7.
#' @seealso \code{\link{tpColorPlot2d}}
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Felsenstein J. (1981). Evolutionary trees from DNA sequences: A maximum likelihood approach. 
#' Journal of Molecular Evolution 17: 368-376.
#'
#' Katoh K. (2013). MAFFT - a multiple sequence alignment program. Available at http://mafft.cbrc.jp/alignment/server.
#'
#' Katoh K, Standley DM. (2013). MAFFT multiple sequence alignment software version 7: 
#' Improvements in performance and usability. Molecular Biology and Evolution 30: 772-780.
#'
#' Minh BQ, Nguyen MAT, Von Haeseler A. (2013). Ultrafast approximation for phylogenetic bootstrap. Molecular Biology and Evolution 30: 1188-1195.
#'
#' Nguyen L-T, Schmidt HA, Von Haeseler A, Minh BQ. (2015). IQ-TREE: a fast and effective stochastic algorithm for estimating maximum likelihood phylogenies. Molecular Biology and Evolution 32: 268-274.
#'
#' Paradis E, Claude J & Strimmer K. (2004). APE: analyses of phylogenetics and evolution in R
#' language. Bioinformatics 20: 289-290.
#' @examples
#' library(ape)
#'
#' data(ligotree)
#' plot.phylo(ligotree, show.node.label=TRUE)
#'
"ligotree"


