#' Read haplotype data, selected by region of interest, from a .haps file
#' (regular ShapeIt output files)
#'
#' @param filename character, path to a \code{.haps} file containing
#' haplotype data which is regular ShapeIt output. See
#' \url{https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#hapsample}
#' @inheritParams read.haplo
#' @return matrix object containing the haplotypes selected by the
#' region of interest
#' @seealso \code{\link{read.haplo}},
#' \code{\link{read.haplo.pedfile}},
#' \code{\link{read.haplo.bedfile}},
#' \code{\link{readMapFile}}
#' @keywords internal
read.haplo.shapeit_haps <- function(filename, map, chr, startpos, endpos){
    # TODO: add extensions
    map <- readMapFile(filename, morgans=FALSE)

    snps2out <- map[which(map[, 1] == chr &
                              map[, 3] > startpos &
                                  map[, 3] < endpos),
                    2]

    write.table(snps2out,
                file = "snps2out.txt",
                quote = FALSE,
		col.names = FALSE,
                row.names = FALSE)

    system(paste0("grep -f snps2out.txt -w ",
                  filename,
                  " >HaploDataSubset.haps"))

    tmpHaplo <- read.table(file = "HaploDataSubset.haps",
                           header = FALSE,
                           stringsAsFactors = FALSE,
                           colClasses = c(rep("NULL", 5),
                               rep("integer",
                                   count.fields(filename)[1] - 5)))

    unlink(c("snps2out.txt", "HaploDataSubset.haps"))

    return(Get.G(tmpHaplo))
}
