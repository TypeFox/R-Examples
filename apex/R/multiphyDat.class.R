
############################
####  CLASSE DEFINITION ####
############################
#'
#' multiphyDat: class for multiple gene data
#'
#' This formal (S4) class is identical to \linkS4class{multidna}, except that DNA sequences are stored using \code{phyDat} objects from the \code{phangorn} package.
#' Sequences are stored as a (possibly named) list, with each element of the list being a separate DNA alignment stored as a \code{phyDat} object.
#' The rows of the separate matrices all correspond to the same individuals, ordered identically.
#'
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @rdname multiphyDat
#'
#' @aliases multiphyDat
#' @aliases multiphyDat-class
#' @aliases listOrNULL
#' @aliases data.frameOrNULL
#'
#' @slot seq a list of phyDat objects; empty slot should be NULL
#' @slot type a character string indicating the type of the sequences stored: "DNA" for DNA sequences, "AA" for amino-acids.
#' @slot labels a vector of labels of individuals
#' @slot n.ind the number of individuals
#' @slot n.seq the total number of sequences (pooling all genes), including gap sequences
#' @slot n.seq.miss the total number of gap-only sequences
#' @slot ind.info a data.frame containing information on the individuals, where individuals are in rows; empty slot should be NULL
#' @slot gene.info a data.frame containing information on the genes, where genes are in rows; empty slot should be NULL
#'
#' @export
#'
#' @import phangorn
#' @import methods
#'
#' @examples
#' data(Laurasiatherian)
#'
#' ## empty object
#' new("multiphyDat")
#'
#' ## simple conversion with nicely ordered output
#' \dontrun{
#' genes <- list(gene1=subset(Laurasiatherian,,1:1600, FALSE),
#'     gene2=subset(Laurasiatherian,,1601:3179, FALSE))
#' x <- new("multiphyDat", genes)
#' x
#' }
#'
#' ## trickier conversion with missing sequences / wrong order
#' genes <- list(gene1=subset(Laurasiatherian,1:40),
#'     gene2=subset(Laurasiatherian,8:47))
#' x <- new("multiphyDat", genes)
#' x
#'
setClass("multiphyDat", representation(seq="listOrNULL", type="character"),
         prototype(seq=NULL, type=character(0)), contains="multiinfo")
