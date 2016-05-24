##' Read pedigree structure from a Plink file
##'
##' Extract the pedigree data from a Plink file (in \code{.ped} or
##' \code{.bed} format). For example, in the case of a \code{.ped}
##' file this simply extracts the first four columns from the file.
##' @param type character, \code{'ped'}, \code{'bed'}: format of the
##' input file (see the \code{filename} parameter) containing the
##' pedigree data (and usually the genomic data as well).
##' @param filename character, path to the input file containing the
##' pedigree data.
##' @return a data frame containing the columns with pedigree data
##' taken from the input file (i.e. family ID, individual ID, father
##' ID, mother ID)
##' @author Lennart C. Karssen
##' @keywords internal
read.pedigree <- function(type, filename) {
    switch(type,
           ped={
               nrcols <- ncol(read.table(file=filename, nrows=2))
               pedigree <- read.table(file=filename,
                                      colClasses=
                                          c(rep("character",
                                                4),
                                            rep("NULL", nrcols - 4)),
                                      stringsAsFactors=FALSE,
                                      header=FALSE
                                      )
           },
           bed={
               filename <- sub(".bed$", ".fam", filename )
               nrcols <- ncol(read.table(file=filename, nrows=2))
               pedigree <- read.table(file=filename,
                                      colClasses=
                                          c(rep("character",
                                                4),
                                            rep("NULL", nrcols - 4)),
                                      stringsAsFactors=FALSE,
                                      header=FALSE
                                      )
           }
           )
    colnames(pedigree) <- c("famID", "indID", "fatherID", "motherID")
    return(pedigree)
}
