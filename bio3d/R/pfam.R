"pfam" <- function(id, alignment='seed', verbose=FALSE) {
  ##alignment <- 'full' ## seed, ncbi, full, metagenomics
  
  oops <- requireNamespace("RCurl", quietly = TRUE)
  if(!oops)
    stop("Please install the RCurl package from CRAN")
  
  cl <- match.call()
  format <- "fasta"
  
  url = paste('http://pfam.sanger.ac.uk/family/', id, '/acc', sep='')
  if(verbose)
    cat("Fetching accession from", url, "\n")
  
  if(!RCurl::url.exists(url)) {
    cat(url, "\n")
    stop("Url does not exist")
  }
  accid <- readLines(url, warn=FALSE)[1]
  
  ## download alignment
  url <- paste('http://pfam.sanger.ac.uk/family/', accid,
               '/alignment/', alignment, '/format?format=', format, sep='')

  if(verbose)
    cat("Fetching alignment from", url, "\n")

  if(!RCurl::url.exists(url)) {
    cat(url, "\n")
    stop("Url does not exist")
  }
  
  tmpfile <- tempfile()
  success <- download.file(url, tmpfile, quiet=!verbose)

  if(success==1)
    stop("Download failed")

  if(verbose)
    cat("Alignment successfully downloaded (", tmpfile, ")\n")
    
  fasta <- read.fasta(tmpfile)
  unlink(tmpfile)
  fasta$call=cl
  return(fasta)
}

