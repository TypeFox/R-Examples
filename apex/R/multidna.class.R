
############################
####  CLASSE DEFINITION ####
############################
setClassUnion("listOrNULL", c("list", "NULL"))

setClassUnion("data.frameOrNULL", c("data.frame", "NULL"))

setClass("multiinfo",
         representation(labels="character", n.ind="integer", n.seq="integer", n.seq.miss="integer",
                        ind.info="data.frameOrNULL", gene.info="data.frameOrNULL", "VIRTUAL"),
         prototype(labels=character(0), n.ind=0L, n.seq=0L, n.seq.miss=0L, ind.info=NULL, gene.info=NULL)
         )

#'
#' multidna: class for multiple gene data
#'
#' This formal (S4) class is used to store multiple DNA alignments.
#' Sequences are stored as a (possibly named) list, with each element of the list being a separate DNA alignment stored as a DNAbin matrix.
#' The rows of the separate matrices all correspond to the same individuals, ordered identically.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @rdname multidna
#'
#' @aliases multidna
#' @aliases multidna-class
#' @aliases listOrNULL
#' @aliases data.frameOrNULL
#'
#' @slot dna a list of DNAbin matrices; empty slot should be NULL
#' @slot labels a vector of labels of individuals
#' @slot n.ind the number of individuals
#' @slot n.seq the total number of sequences (pooling all genes), including gap sequences
#' @slot n.seq.miss the total number of gap-only sequences
#' @slot ind.info a data.frame containing information on the individuals, where individuals are in rows; empty slot should be NULL
#' @slot gene.info a data.frame containing information on the genes, where genes are in rows; empty slot should be NULL
#'
#' @export
#'
#' @import ape
#' @import methods
#'
#' @examples
#'
#' ## empty object
#' new("multidna")
#'
#' ## simple conversion with nicely ordered output
#' data(woodmouse)
#' genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
#' x <- new("multidna", genes)
#' x
#' image(woodmouse)
#' image(x@@dna[[1]])
#' image(x@@dna[[2]])
#'
#' ## trickier conversion with missing sequences / wrong order
#' genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[c(5:1,14:15),501:965])
#' x <- new("multidna", genes)
#' x
#' image(x@@dna[[1]])
#' image(x@@dna[[2]])
#'
setClass("multidna", representation(dna="listOrNULL"), prototype(dna=NULL), contains="multiinfo")

