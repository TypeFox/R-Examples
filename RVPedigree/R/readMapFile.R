#' Read file with information about SNPs chromosome and position, for example
#' from regular PLINK .map OR .bim file
#'
#' These files have (at least) the following columns (separated by
#' white space):
#' \itemize{
#' \item column 1: chromosome
#' \item column 2: variant name
#' \item column 3: genetic distance in morgans (optional, see the
#' \code{morgans} option)
#' \item column 4: base-pair position (bp units)
#' }
#' @param filename character, path to input file containing genomic
#' map data, e.g. a Plink
#' \href{http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map}{\code{.map}},
#' \href{http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#tr}{\code{.bim}}
#' or \href{http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed}{\code{.tped}} file.
#' @param morgans logical, indicate whether input file contains column with
#' genetic distance between SNPs (in which case the bp positions are
#' in the fourth column)
#' @return matrix object with columns contains chromosome, rsID
#' and position in bp
#' @seealso \code{\link{read.haplo}},
#' \code{\link{read.haplo.pedfile}},
#' \code{\link{read.haplo.bedfile}},
#' \code{\link{read.haplo.shapeit_haps}}
#' @author Sodbo Sharapov
#' @export
readMapFile <- function(filename = "NULL", morgans = TRUE){
    # TODO: change count.fields() to length(readLine(file,n=1))
    # TODO: add extensions
    if (!morgans) {
        mapInfo <- read.table(file = filename,
                              colClasses =
                                  c("character",
                                    "character",
                                    "numeric",
                                    rep("NULL",
                                        count.fields(filename)[1] - 3)),
                              stringsAsFactors = FALSE,
                              header = FALSE)
    } else {
        mapInfo <- read.table(file = filename,
                              colClasses =
                                  c("character",
                                    "character",
                                    "NULL",
                                    "numeric",
                                    rep("NULL",
                                        count.fields(filename)[1] - 4)),
                              stringsAsFactors = FALSE,
                              header = FALSE)
    }
    return(mapInfo)
}

# tmp1 <- readMapFile(filename = "data.map",morgans = TRUE)
# head(tmp1)
# tmp2 <- readMapFile(filename = "dataB.bim",morgans = TRUE)
# head(tmp2)
# tmp3 <- readMapFile(filename = "dataT.tped",morgans = TRUE)
# head(tmp3)
# tmp4 <- readMapFile(filename ="data_nullmorgans.map",morgans = FALSE)
# head(tmp4)
