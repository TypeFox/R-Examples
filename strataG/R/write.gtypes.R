#' @title Write \code{gtypes}
#' @description Write a \linkS4class{gtypes} object to file(s).
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param label label for filename(s). Default is the gtypes description 
#'   if present.
#' @param folder folder where file(s) should be written to. If \code{NULL}, 
#'   files are written to current working directory.
#' @param as.frequency logical indicating if haploid data should be output 
#'   as frequency tables.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom utils write.csv
#' @export
#' 
write.gtypes <- function(g, label = NULL, folder = NULL, as.frequency = FALSE) {
  desc <- description(g)
  label <- if(!is.null(label)) {
    label 
  } else if(!is.null(desc)) {
    desc 
  } else "strataG.gtypes"
  label <- gsub("[[:punct:]]", ".", label)
  out.files <- paste(label, ".csv", sep = "")
  if(!is.null(folder)) out.files <- file.path(folder, out.files)
  
  g.mat <- if(ploidy(g) == 1 & as.frequency) {
    x <- as.frequency(g) 
    x <- data.frame(haplotype = rownames(x), cbind(x))
    rownames(x) <- NULL
    x
  } else as.matrix(g)
  g.mat <- cbind(id = rownames(g.mat), strata = strata(g), g.mat)
  write.csv(g.mat, file = out.files, row.names = FALSE)
  
  if(!is.null(sequences(g))) {
    for(x in locNames(g)) {
      fname <- paste(label, x, "fasta", sep = ".")
      if(!is.null(folder)) fname <- file.path(folder, fname)
      write.dna(sequences(g, x), file = fname, format = "fasta", nbcol = -1, 
                colsep = "", indent = 0, blocksep = 0)
      out.files <- c(out.files, fname)
    }
  }
  
  invisible(out.files)
}