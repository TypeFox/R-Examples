# some auxiliary functions from package 'seqminer' (c) 2015

#' Check if the input is url e.g. http:// or ftp://
#' @param fileName character vector
#' @keywords internal 

isURL <- function(fileName) {
  if (grepl(pattern = "^http://", fileName) ||
      grepl(pattern = "^ftp://", fileName) ) {
    return(TRUE)
  }
  return(FALSE)
}

local.file.exists <- function(fileName) {
  if (isURL(fileName)) {
    return(TRUE)
  }
  return(file.exists(fileName))
}

#' Check input file has tabix index
#'
#' @param fileName character vector, representing file names an input file name
#' @return TRUE if an index file with ".tbi" or "bci" exists, and created later than VCF/BCF file
#' @keywords internal

hasIndex <- function(fileName) {
  if (isURL(fileName)) {
    return(TRUE)
  }
  endsWith <- function(i, pattern = NULL) {
    if (is.null(pattern)) {
      return(FALSE)
    }
    if (length(fileName) != 1) {
      stop("hasIndex() only take vector of length 1.")
      return (FALSE)
    }
    p <- paste(pattern, "$", sep = "")
    if (length(grep(x=i, pattern=p)) != 1) {
      return (FALSE)
    }
    return (TRUE)
  }
  ## handle VCF
  fIndex <- NULL
  if (endsWith(fileName, ".vcf")) {
    stop(gettextf("Input file is not bgzip compressed. Please use: [ bgzip %s ]", fileName))
    return (FALSE)
  }
  if (endsWith(fileName, ".vcf.gz")) {
    fIndex <- paste(fileName, ".tbi", sep = "")
    if (!file.exists(fIndex)) {
      stop(gettextf("Cannot find index file (you can create it using: [ tabix -p vcf %s ])" , fileName))
      return (FALSE)
    }
  }

  ## handle BCF
  if (endsWith(fileName, ".bcf")) {
    stop(gettextf("Input file is not bgzip compressed. Please use: [ bgzip %s ]", fileName))
    return (FALSE)
  }
  if (endsWith(fileName, ".bcf.gz")) {
    fIndex <- paste(fileName, ".bci", sep = "")
    if (!file.exists(fIndex)) {
      stop(gettextf("Cannot find index file (you can create it using: [ bcftools index %s ])" , fileName))
      return (FALSE)
    }
  }

  ## Check index file
  if (is.null(fIndex)) {
    stop(gettextf("Cannot guess a valid tabix index, did your file have suffix: .vcf, .vcf.gz, .bcf, or .bcf.gz ?"))
    return (FALSE)
  }
  if (!file.exists(fIndex)) {
    stop(gettextf("The index file [ %s ] does not exists.", fIndex))
    return (FALSE)
  }

  ## Check index file timestamp
  time.data <- file.info(fileName)$mtime
  time.index <- file.info(fIndex)$mtime
  if (time.data > time.index) {
    stop(gettextf("Index file [ %s ] is older than data file [ %s ].", fIndex, fileName))
    return(FALSE)
  }
  return (TRUE)
}
