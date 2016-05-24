##' wc
##' word count of ASCII files
##' 
##' wc uses the system command wc
##' 
##' @param file the file name or pattern
##' @param flags the parameters to count, character vector with the long form
##'   of the parameters
##' @return data.frame with the counts and file names, or \code{NULL} if wc is
##'   not available
##' @export
##' @author C. Beleites
wc <- function (file, flags = c("lines", "words", "bytes")){
  if (length (system ("wc --help", intern = TRUE)) == 0)
    return (NULL)

  wc <- paste ("wc", paste ("--", flags, sep = "", collapse = ", "), file)
  wc <- read.table(pipe (wc))
  colnames (wc) <- c(flags, "file")
  wc
}

