#'
#' Convert multidna into genind
#'
#' The functions \code{multidna2genind} and \code{multiphyDat2genind} concatenates separate DNA alignments, and then extracts SNPs of the resulting alignment into a \linkS4class{genind} object.
#'
#' @docType methods
#' @rdname multidna2genind
#'
#' @param x a \linkS4class{multidna} or \linkS4class{multiphyDat} object.
#' @param genes an optional vector indicating the genes to retain for the concatenation; any way to subset the list in x@@dna or x@@seq is acceptable; by default, all genes are used.
#' @param mlst if \code{TRUE}, each gene will result in a single locus in the genind object. (Default to \code{FALSE})
#' @param gapIsNA if \code{TRUE} and \code{mlst = TRUE}, sequences that consist entirely of gaps will be considered as NAs. (Default to \code{FALSE})
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}, Zhian N. Kamvar, Klaus Schliep
#'
#' @rdname multidna2genind
#'
#' @aliases multidna2genind
#' @aliases multiphyDat2genind
#'
#' @seealso
#' \itemize{
#' \item concatenate
#' \item \code{\link{DNAbin2genind}} to convert single DNAbin objects.
#' }
#'
#' @export
#'
#' @importFrom adegenet DNAbin2genind df2genind
#'
#' @return a \linkS4class{genind} object
#'
#' @examples
#'
#' ## simple conversion with nicely ordered output
#' data(woodmouse)
#' genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
#' x <- new("multidna", genes)
#' x
#' y <- multidna2multiphyDat(x)
#' y
#' z1 <- multidna2genind(x)
#' z1
#' z2 <- multiphyDat2genind(y)
#' all.equal(z1, z2)
#'
multidna2genind <- function(x, genes=TRUE, mlst=FALSE, gapIsNA=FALSE){
  if (!mlst){
    return(DNAbin2genind(concatenate(x, genes=genes)))
  }
  xlist  <- lapply(x@dna, function(i) apply(as.character(i), 1, paste, collapse = ""))
  xdf    <- data.frame(xlist)
  xlevs  <- lapply(xdf, levels)
  if (gapIsNA){
    xgap <- find_gap_sequence(xlevs)
    for (i in names(xgap)){
      the_gap <- xgap[[i]]
      levels(xdf[[i]])[the_gap] <- NA
      if (length(the_gap) > 0){
        xlevs[[i]] <- xlevs[[i]][-the_gap]
      }
    }
  }
  xdfnum <- data.frame(lapply(xdf, as.numeric))
  xgid   <- df2genind(xdfnum, ploidy = 1, ind.names = x@labels)
  names(xlevs)   <- names(xgid@all.names)
  xgid@all.names <- xlevs
  xgid@other$ind.info <- x@ind.info
  xgid@other$gene.info <- x@gene.info
  return(xgid)
}


#' @rdname multidna2genind
#' @export
#'
multiphyDat2genind <- function(x, genes=TRUE, mlst=FALSE, gapIsNA=FALSE){
    return(multidna2genind(multiphyDat2multidna(x), genes=genes, mlst=mlst, gapIsNA=gapIsNA))
}





#'
#' Conversions between multidna and multiphyDat
#'
#' The functions \code{multidna2multiphyDat} and \code{multiphyDat2multidna} are used to convert data between \linkS4class{multidna} and  \linkS4class{multiphyDat} classes.
#'
#' @docType methods
#'
#' @param x a \linkS4class{multidna} or \linkS4class{multiphyDat} object.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}, Zhian N. Kamvar, Klaus Schliep
#'
#' @rdname multidna2multiphyDat
#'
#' @aliases multidna2multiphyDat
#' @aliases multiphyDat2multidna
#'
#' @seealso
#' \itemize{
#' \item concatenate
#' \item \code{\link{DNAbin2genind}} to convert single DNAbin objects.
#' }
#'
#' @export
#'
#' @examples
#'
#' ## simple conversion with nicely ordered output
#' data(woodmouse)
#' genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
#' x <- new("multidna", genes)
#' x
#'
#' ## conversion multidna -> multiphyDat
#' y <- multidna2multiphyDat(x)
#' y
#'
#' ## check round trip
#' identical(x, multiphyDat2multidna(y))
#'
#' @rdname multidna2multiphyDat
#'
#' @export
#'
multidna2multiphyDat <- function(x){
    tmp <- lapply(x@dna, phyDat)
    new("multiphyDat",tmp, type="DNA")
}


#' @rdname multidna2multiphyDat
#' @export
multiphyDat2multidna <- function(x){
    if(!tolower(x@type) %in% c("dna","rna")) stop("Conversion to multidna only possible with DNA sequences ('@type' should be 'DNA')")
    tmp <- lapply(x@seq, as.character)
    new("multidna",tmp)
}





#'
#' Convert from multidna into alignment (seqinr)
#'
#' The functions \code{multidna2alignment} and \code{multiphyDat2alignment} concatenates separate sequences and return an alignment object of the seqinr package.
#'
#' @docType methods
#'
#' @param x a \linkS4class{multidna} or \linkS4class{multiphyDat} object.
#' @param genes an optional vector indicating the genes to retain for the concatenation; any way to subset the list in x@@dna or x@@seq is acceptable; by default, all genes are used.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}, Zhian N. Kamvar, Klaus Schliep
#'
#' @rdname multidna2alignment
#'
#' @aliases multidna2alignment
#' @aliases multiphyDat2alignment
#'
#' @seealso
#' \itemize{
#' \item concatenate
#' \item \code{\link{as.alignment}} to convert single DNAbin objects.
#' }
#'
#' @export
#'
#' @return a alignment object
#'
#' @examples
#' \dontrun{
#' ## simple conversion with nicely ordered output
#' data(woodmouse)
#' genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
#' x <- new("multidna", genes)
#' x
#' y <- multidna2alignment(x)
#' y
#' x2 <- multidna2multiphyDat(x)
#' z <- multiphyDat2alignment(x2)
#' }
#' @export
#'
multidna2alignment <- function(x, genes=TRUE){
  return(as.alignment(concatenate(x, genes)))
}


#' @rdname multidna2alignment
#' @export
#'
multiphyDat2alignment <- function(x, genes=TRUE){
  return(phyDat2alignment(concatenate(x, genes)))
}





##
## internal functions
##
find_gap_sequence <- function(x){
  wheregaps <- lapply(x, function(i) which(grepl("^\\-+?$", i)))
  return(wheregaps)
}

