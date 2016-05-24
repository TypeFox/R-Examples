#' Sample from the rows of a (possibly large) text file (NOT IMPLEMENTED)
#'
#' Read a large text file in batches, keeping the rows to be included in the sample. Return a data.frame whose rows are the sample.
#' 
#' @param n numeric, size of sample to be taken.
#' @param out.of numeric, number of rows in the data set not including the header.
#' @param filename character, name of the file containing the data. This must be a tab-delimited file with a header row formatted per the default options for \code{\link{read.delim}}.
#' @return data.frame, size n random subset of the rows of filename
#' @seealso \code{\link{ScalablePCA}}, \code{\link{GetSampleFromDataFrame}}, \code{\link{GetSampleFromDb}}
#' @references
#' \url{https://github.com/shaptonstahl/}
#' @author Stephen R. Haptonstahl \email{srh@@haptonstahl.org}
#' @examples
#' \dontrun{x <- dils:::GetSampleFromFile(10, 150, "folder/containing/data.txt")}
GetSampleFromFile <- function(n,
                              out.of,
                              filename) {  
  # initialize output
  out <- read.delim(filename, nrows=5)
  out <- out[numeric(0),]
  
  rows.to.get <- sort(sample(1:out.of, n))
  
  con <- file(filename, "r")
  garbage <- readLines(con, n=1); rm(garbage)  # read and toss header
  n.rows.to.read <- 1e4
    
  while( length(rows.to.get) > 0 ) {
    input.text <- readLines(con, n=n.rows.to.read)  # no check for empty input.text
    
    use.indices <- rows.to.get[rows.to.get <= n.rows.to.read]
    
    if( length(use.indices > 0) ) {
      use.text <- input.text[use.indices]
      new.df <- read.delim(header=FALSE, text=use.text)
      out <- rbind(out, new.df)
    }
    rows.to.get <- rows.to.get - n.rows.to.read
    rows.to.get <- rows.to.get[rows.to.get > 0]
  }
  close(con)
  
  return(out)
}