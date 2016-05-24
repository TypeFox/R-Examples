#' Convert PopGenome Data into Coala's Format
#'
#' Using this function, you can convert genetic data imported with the package
#' \pkg{PopGenome} into coala's segsites format. See \code{\link{as.segsites}}
#' for general information on converting genetic data for coala.
#'
#' This function imports all loci from the \code{GENOME} object that have at
#' least one valid site (\code{data@n.valid.sites}). The number of valid sites
#' is used as length of a locus.
#'
#' @param data The \code{GENOME} data from \pkg{PopGenome}.
#' @param only_synonymous Only import synonymous SNPs if set to \code{TRUE}.
#'   This requires that \pkg{PopGenome} knows where coding regions are., e.g.
#'   by using gff files.
#' @param ... Ignored.
#' @seealso An example and additional instructions are provided in the
#'          `coala-data-import` vignette
#' @export
as.segsites.GENOME <- function(data, only_synonymous = FALSE, ...) {
  require_package("PopGenome")

  lapply(which(data@n.valid.sites > 0), function(locus_number) {
    bam <- PopGenome::get.biallelic.matrix(data, locus_number)

    if (is.null(bam)) {
      warning("Locus ", locus_number, " is NULL")
      return(create_empty_segsites(length(unlist(data@populations)) +
                                     length(data@outgroup)))
    }

    pops <- which(row.names(bam) %in% unlist(data@populations))
    outgroup <- which(row.names(bam) %in% data@outgroup)
    stopifnot(!is.null(pops))
    stopifnot(!is.null(outgroup))

    # Select relevant data
    if (only_synonymous) {
      syn <- data@region.data@synonymous[[locus_number]]
      if (is.null(syn)) syn <- numeric(0)
      syn[is.na(syn)] <- FALSE
      seg_sites <- bam[c(pops, outgroup), syn, drop = FALSE]
    } else {
      seg_sites <- bam[c(pops, outgroup), , drop = FALSE]
    }
    stopifnot(!is.null(seg_sites))

    # Add positions attribute
    pos <- as.numeric(colnames(bam)) / data@n.sites[[locus_number]]
    create_segsites(seg_sites, pos)
  })
}


#' @importFrom utils capture.output
create_popgenome_test_data <- function() {
  require_package("PopGenome")

  capture.output({
    fasta <- system.file("example_fasta_files", package = "coala")
    data_pg <- PopGenome::readData(fasta, progress_bar_switch = FALSE)
    data_pg <- PopGenome::set.outgroup(data_pg, c("Individual_Out-1",
                                                  "Individual_Out-2"))
    data_pg <- PopGenome::set.populations(data_pg,
                                          list(paste0("Individual_1-", 1:5),
                                               paste0("Individual_2-", 1:5)))
  })
  data_pg
}
