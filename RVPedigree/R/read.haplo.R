#' Read haplotype data, selected by region of interest, from PLINK formatted
#' files or ShapeIt output files
#'
#' @param type character, \code{'ped'}, \code{'bed'} (default) or
#' \code{'shapeit-haps'} format of input file containing haplotype
#' data
#' @param filename character, path to input file containing haplotype data
#' @param map object, data.frame contains 3 columns: rsID, chromosome,
#' position in bp as output by e.g. \code{\link{readMapFile}}.
#' @param chr character, chromosome number (basically from 1 to 22 as used by
#' \href{http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped}{Plink}),
#' on which the region of interest is located
#' @param startpos numeric, start position (in bp, base pairs) of the
#'     region of interest (default: 0)
#' @param endpos numeric, end position (in bp, base pairs) of the
#'     region of interest (default: 0)
#' @return matrix object containing the haplotypes selected by the
#' region of interest
#' @seealso \code{\link{read.haplo.pedfile}},
#' \code{\link{read.haplo.bedfile}},
#' \code{\link{read.haplo.shapeit_haps}},
#' \code{\link{readMapFile}}
#' @keywords internal
read.haplo <- function(type="bed",
                       filename,
                       map,
                       chr=0,
                       startpos=0,
                       endpos=0){
    check_positions(startpos=startpos, endpos=endpos)
    ## TODO: add extensions
    if (type == 'bed') {
        return (read.haplo.bedfile(filename = filename,
                                   map = map,
                                   chr = chr,
                                   startpos = startpos,
                                   endpos = endpos))
    }

    if (type == 'ped') {
        return (read.haplo.pedfile(filename = filename,
                                   map = map,
                                   chr = chr,
                                   startpos = startpos,
                                   endpos = endpos))
    }

    if (type == '') {
        return (read.haplo.shapeit_haps(filename = filename,
                                        map = map,
                                        chr = chr,
                                        startpos = startpos,
                                        endpos = endpos))
    }
}

# source("../../../shapeit/R/read.haplo.pedfile.R")
# source("../../../shapeit/R/read.haplo.bedfile.R")
# source("../../../shapeit/R/readMapFile.R")

# plinkMap <- readMapFile(filename = "data.map",morgans = TRUE)

# PED <- read.haplo(type = "ped", filename = "data.ped", map = plinkMap,
#		chr = 1, startpos = 9000, endpos = 25000)

# BED <- read.haplo(type = "bed", filename = "dataB", map = plinkMap,
#		chr = 1, startpos = 9000, endpos = 25000)
