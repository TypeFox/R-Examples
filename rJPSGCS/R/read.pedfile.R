# Adapted from snpMatrix read.snps.pedfile.R

read.pedfile.info <- function(file) {
  # do exactly 3 columns, even if there are less or more
  result <- read.table(file, col.names=c("snp.names", "position", "chromosome"),
                       as.is=TRUE, fill=TRUE, flush=TRUE)
  rownames(result) <- result$snp.names
  result
}

read.pedfile.map <- function(file) {
  # return exactly 4 columns, even if there are less or more
  result <- read.table(file, col.names=c("chromosome", "snp.names", "genetic.distance", "position"),
                       as.is=TRUE, fill=TRUE, flush=TRUE)
  rownames(result) <- result$snp.names
  result[,c('snp.names','position', 'chromosome')] # drop genetic distance and reorder cols for consistency
}

read.pedfile <- function(file, snp.names=NULL, assign=NULL, missing=NULL, 
  X=FALSE, sep=".") {
  low.mem=FALSE
  # If there is no input names, try to see if we can find the accompanying info file,
  # and if possible, load it for the snp names
  join.info <- FALSE
  if (is.null(snp.names)) {
    map.file  <- paste(file, "map", sep=".")
    info.file <- paste(file, "info", sep=".")
    if(length(grep("\\.ped$", file))) {
      map.file  <- sub('\\.ped$', '.map',  file)
      info.file <- sub('\\.ped$', '.info', file)
    }
    if (!(file.access(map.file,mode=4))) {
      #cat("Found accompanying map file, reading it first\n")
      snp.info <- read.pedfile.map(map.file)
      snp.names <- rownames(snp.info)
      join.info <- TRUE    
    } else
    if (!(file.access(info.file,mode=4))) {
      # file.access() return 0 for success
      #cat("Found accompanying info file, reading it first\n")
      snp.info <- read.pedfile.info(info.file)
      snp.names <- rownames(snp.info)
      join.info <- TRUE
    }
  }
  if (!is.null(snp.names)) {
    snp.names <- as.character(snp.names)
  }
  if (!is.null(missing)) {
    missing <- as.character(missing)
  }
  
  # if file starts with http:// or ftp://, download it and replace the input
  # with the downloaded file
  if (length(grep("^(ftp|http|file)://", file)) > 0) { 
  # mode = b is needed for windows
    saved.file <- tempfile()
    status <- download.file(file, destfile=saved.file, mode="wb")    
  # download.file() supposedly should throw error already, or return 0/1
    if ((status != 0) && (status != 1)) {
      stop("Download has gone wrong\n");
    }
    file <- saved.file
  }
  
  result <- .Call("read_pedfile", file, snp.names, missing,
                  as.logical(X), as.character(sep), PACKAGE="rJPSGCS")
  if (join.info) {
    snp.info[['assignment']] <- as.factor(result$snp.support)
    result$snp.support <- snp.info
  } else {
    result$snp.support <- as.factor(result$snp.support)
  }

  result
}
