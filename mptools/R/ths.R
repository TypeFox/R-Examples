#' Extract THS from a folder of RAMAS Metapop .ptc files
#' 
#' Extracts Total Habitat Suitability from RAMAS .ptc file.
#' 
#' @param ptc A character string containing the path to a directory containing
#'    the .ptc files of interest. This should contain no .ptc other than those
#'    to be included in the extraction. Order will be taken from the numeric 
#'    component of filenames.
#' @return A numeric vector with length equal to the number of .ptc files in the 
#'    directory specified by \code{ptc}.
#' @seealso \code{\link{results}}
#' @note This has been tested for RAMAS version 5, and may produce unexpected
#'   results for other versions.
#' @importFrom utils tail
#' @export
ths <- function(ptc) {
  ptcs <- list.files(ptc, pattern='\\.ptc$', ignore.case=TRUE, full.names=TRUE)
  ptcs <- ptcs[order(as.numeric(gsub('\\D', '', basename(ptcs))))]
  get.ths <- function(ptc) {
    txt <- readLines(ptc)[-1]
    if(any(grepl('^Results:', txt))) {
      txt <- txt[-(1:grep('^Results:', txt))]
      txt <- txt[-(1:(
        utils::tail(which(sapply(gregexpr(',', txt), length) == 25), 1) + 2))]
      txt <- txt[1:(which(sapply(gregexpr(' ', txt), length) < 6)[1] - 1)]
      hs <- apply(do.call(rbind, strsplit(txt, ' ')), 2, as.numeric)
    } else {
      hs <- matrix(0, ncol=7)
    }
    if (!is.matrix(hs)) hs[1] else sum(hs[, 1])
  }
  message('Calculating total habitat suitability from ptc files in:\n', ptc)
  sapply(ptcs, get.ths, USE.NAMES=FALSE)
}
