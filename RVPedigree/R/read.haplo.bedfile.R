#' Read haplotype data, selected by region of interest, from BED files
#' (regular PLINK binary files)
#'
#' @param filename character, path to BED file containing haplotype
#' data. Attention! filename should contain path to \code{<file>.bed}
#' with full file name. For example \code{../mydata/inputplink.bed}.
#' @inheritParams read.haplo
#' @return matrix object containing the haplotypes selected by the
#' region of interest
#' @seealso \code{\link{read.haplo}},
#' \code{\link{read.haplo.pedfile}},
#' \code{\link{read.haplo.shapeit_haps}},
#' \code{\link{readMapFile}}
#' @import snpStats
#' @keywords internal
read.haplo.bedfile <- function(filename = "NULL",
                               map,
                               chr = "NULL",
                               startpos = "NULL",
                               endpos = "NULL"){
    # TODO: add extensions
    snps2out <- map[which(map[, 1] == chr &
                              map[, 3] > startpos &
                                  map[, 3] < endpos),
                    2]

    plink.input <- snpStats::read.plink(bed = paste0(filename, ".bed"),
                                        bim = paste0(filename, ".bim"),
                                        fam = paste0(filename, ".fam"),
                                        select.snps = snps2out)

    # TODO: added rownames(GenotyMatrix) <- idNames may requied
    GenotypeMatrix <- methods::as(plink.input$genotypes, "numeric")
    return(GenotypeMatrix)
}

# tmp = read.haplo.bedfile(filename = "dataB",map = tmp1, chr = 1,
#	startpos = 9000, endpos = 25000)
