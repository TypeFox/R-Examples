
#' Read multiple DNA alignments
#'
#' These functions read multiple DNA alignments and store the output in a \linkS4class{multidna} object.
#' They are relying on ape's original functions \code{\link[ape]{read.dna}} and \code{\link[ape]{read.FASTA}}.
#'
#' @rdname readfiles
#' @aliases read.multidna
#' @aliases read.multiFASTA
#' @aliases read.multiphyDat
#'
#'
#' @param files a vector of characters indicating the paths to the files to read from.
#' @param add.gaps a logical indicating if gap-only sequences should be added wherever sequences are missing; defaults to TRUE.
#' @param ... further arguments to be passed to the functions \code{\link[ape]{read.dna}} and \code{\link[ape]{read.FASTA}}.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#'
#' @seealso
#' \itemize{
#' \item \code{\link[ape]{read.dna}}
#' \item  \code{\link[ape]{read.FASTA}}
#' \item \code{\link[phangorn]{read.phyDat}}
#' }
#'
#' @export
#'
#' @examples
#' ## get path to the files
#' files <- dir(system.file(package="apex"),patter="patr", full=TRUE)
#' files
#'
#' ## read files
#' x <- read.multiFASTA(files)
#' x
#' plot(x)
#'
#' y <- read.multiphyDat(files, format="fasta")
#' y
#'
read.multidna <- function(files, add.gaps=TRUE, ...){
    gene.names <- gsub(".fasta","", basename(files))
    dna <- lapply(files, read.dna, ...)
    names(dna) <- gene.names
    out <- new("multidna", dna=dna, add.gaps=add.gaps)
    return(out)
}


#'
#' @rdname readfiles
#' @export
read.multiFASTA <- function(files, add.gaps=TRUE){
    gene.names <- gsub(".fasta","",basename(files))
    dna <- lapply(files, read.FASTA)
    names(dna) <- gene.names
    out <- new("multidna", dna=dna, add.gaps=add.gaps)
    return(out)
}


#'
#' @rdname readfiles
#' @export
read.multiphyDat <- function(files, add.gaps=TRUE, ...){
  gene.names <- gsub(".fasta","",basename(files))
  seq <- lapply(files, read.phyDat, ...)
  names(seq) <- gene.names
  out <- new("multiphyDat", seq=seq, add.gaps=add.gaps)
  return(out)
}
