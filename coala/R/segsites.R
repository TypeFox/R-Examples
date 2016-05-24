#' Segregating Sites
#'
#' These functions create and modify segregating sites
#' objects, which are one of the basic intermediary statistics that is
#' calculated in coala. Segregating sites consist primarily of a SNP matrix that
#' contains all SNPs for one locus, with some additional information attached.
#' The parts of the S3 class are detailed below.
#'
#' A segregating sites object contains all SNPs for one genetic locus. Each
#' object consists of three parts: A SNP matrix, a vector of SNP positons and
#' a vector that states which transcript a SNP belong to, if the locus
#' consists of multiple transscripts ('locus trio').
#'
#' \itemize{
#'   \item{In the \strong{SNP} matrix, each row represents a haplotype and each
#'         column represents a SNP. An entry is either \code{1} if the
#'         haplotype carries the derived allele for the SNP, or \code{0} if it
#'         carries the ancestral one.}
#'   \item{In the \strong{positions} vector, each entry gives the relative
#'         position of SNP in the corresponding column of the SNP matrix.}
#'   \item{The \strong{trio_locus} vector contains the trio locus each SNP
#'         belongs to. Entry of \code{-1},\code{0}, \code{1} represent the
#'         left, middle, and right locus, respectively. For normal loci,
#'         this just consists of \code{0}'s}
#' }
#'
#' @name segsites
#' @author Paul Staab
#' @seealso For converting biological data to segsites: \code{\link{as.segsites}}
#' @examples
#' snps <- matrix(c(1, 0, 0,
#'                  1, 1, 0,
#'                  0, 0, 1), 3, 3, TRUE)
#' pos <- c(.1, .3, .5)
#' segsites <- create_segsites(snps, pos)
#' print(segsites)
#' get_snps(segsites)
#' get_positions(segsites)
#'
#' # When subsetting individuals, sites that are not
#' # segregating in these are automatically removed:
#' segsites[1:2, 1:3]
NULL

#' @export
print.segsites <- function(x, ...) {
  snps <- get_snps(x)
  colnames(snps) <- format(get_positions(x), scientific = FALSE)
  print(snps)
}


#' @export
"[.segsites" <- function(x, chrs, snps, drop = FALSE) {
  create_segsites(snps = get_snps(x)[chrs, select = snps, drop = FALSE],
                  positions = get_positions(x)[snps],
                  trio_locus = get_trio_locus(x)[snps])
}


#' @export
as.matrix.segsites <- function(x, ...) get_snps(x)

#' @export
dim.segsites <- function(x) dim(get_snps(x))



#' @describeIn segsites Checks whether an object is a segsites object.
#' @export
is_segsites <- function(segsites) inherits(segsites, "segsites")


create_test_segsites <- function() {
  create_segsites((matrix(c(1, 1, 0, 1, 1,
                            0, 0, 0, 1, 0,
                            0, 1, 1, 0, 1), 3, 5, byrow = TRUE)),
                  c(.1, .2, .5, .7, .75))
}


create_empty_segsites <- function(n_ind = 0) {
  create_segsites(matrix(0, n_ind, 0), numeric(0))
}


conv_to_ms_output <- function(segsites) {
  c(paste("segsites:", ncol(segsites)),
    paste("positions:", paste(format(get_positions(segsites),
                                     scientific = FALSE),
                              collapse = " ")),
    apply(segsites, 1, paste, collapse = ""))
}



#' @describeIn segsites Combines three segregating sites to a locus trio
#'
#' @param left The segregating sites from the left locus
#' @param middle The segregating sites from the middle locus
#' @param right The segregating sites from the right locus
create_locus_trio <- function(left, middle, right) {
  assert_that(is.list(left) && is.list(middle) && is.list(right))
  assert_that(length(left) == length(middle))
  assert_that(length(left) == length(right))

  lapply(seq(along = left), function(locus) {
    create_segsites(cbind(get_snps(left[[locus]]),
                          get_snps(middle[[locus]]),
                          get_snps(right[[locus]])),
                    c(get_positions(left[[locus]]),
                      get_positions(middle[[locus]]),
                      get_positions(right[[locus]])),
                    c(rep(-1, ncol(left[[locus]])),
                      rep( 0, ncol(middle[[locus]])),
                      rep( 1, ncol(right[[locus]]))),
                    check = FALSE)
  })
}


#' Convert genetic data to coala's internal format
#'
#' This function can be used to convert the genomic data formats used in other
#' packages to calculate coala's segregating sites object. This is useful for
#' calculating the summary statistics of biological data using the
#' \code{\link{calc_sumstats_from_data}} function.
#'
#' Currently, only the package \pkg{PopGenome} is supported, see
#' \code{\link{as.segsites.GENOME}} for details.
#'
#' @param data The data object that is converted.
#' @param ... Additional arguments specific for the used package.
#'
#' @return A list of segregating sites objects.
#' @seealso Further instructions are provided in the `coala-data-import` vignette
#' @seealso For information about segsites: \code{\link{segsites}}
#' @export
as.segsites <- function(data, ...) UseMethod("as.segsites")

#' @export
as.segsites.default <- function(data, ...) {
  stop("Unknown data format")
}

#' @export
as.segsites.list <- function(data, ...) {
  assert_that(all(vapply(data, is_segsites, logical(1))))
  data
}
