#' Read haplotype data, selected by region of interest, from PED files
#' (regular PLINK untransposed text files)
#'
#' @param filename character, path to PED file containing haplotype
#' data. Attention! filename should contain path to \code{<file>.ped}
#' with full file name. For example \code{../mydata/inputplink.ped}
#' @inheritParams read.haplo
#' @return matrix object containing the haplotypes selected by the
#' region of interest
#' @seealso \code{\link{read.haplo}},
#' \code{\link{read.haplo.bedfile}},
#' \code{\link{read.haplo.shapeit_haps}},
#' \code{\link{readMapFile}}
#' @import snpStats
#' @keywords internal
read.haplo.pedfile <- function(filename = "NULL",
                               map,
                               chr = "NULL",
                               startpos = "NULL",
                               endpos = "NULL"){
    # TODO: probably we want to automatically append ".ped" to filename
    # TODO: add extensions
    snps2out <- map[which(map[, 1] == chr &
                              map[, 3] > startpos &
                                  map[, 3] < endpos),
                    2]

    plink.input <- snpStats::read.pedfile(file = filename,
                                          snps = snps2out)

    GenotypeMatrix <- methods::as(plink.input$genotypes, "numeric")
    # TODO: added 'rownames(GenotypeMatrix) <- idNames' may required
    return(GenotypeMatrix)
}

# tmp = read.haplo.pedfile(filename = "data.ped",map = tmp1,chr = 1,
#	startpos = 9000, endpos = 25000)
