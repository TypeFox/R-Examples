## Compress and replace file on computer
xzcompress <- function (filename, destname = sprintf("%s.xz", filename), temporary = FALSE, 
          skip = FALSE, overwrite = FALSE, remove = TRUE, BFR.SIZE = 1e+07, compression = 6,
          ...) 
{
  if (!file.exists(filename)) {
    stop("No such file: ", filename)
  }
  if (temporary) {
    destname <- file.path(tempdir(), basename(destname))
  }
  attr(destname, "temporary") <- temporary
  if (filename == destname) 
    stop(sprintf("Argument 'filename' and 'destname' are identical: %s", 
                 filename))
  if (file.exists(destname)) {
    if (skip) {
      return(destname)
    }
    else if (!overwrite) {
      stop(sprintf("File already exists: %s", destname))
    }
  }
  destpath <- dirname(destname)
  if (!file.info(destpath)$isdir) 
    dir.create(destpath)
  inn <- file(filename, open = "rb")
  on.exit(if (!is.null(inn)) close(inn))
  outComplete <- FALSE
  out <- xzfile(destname, open = "wb", compression = compression, ...)
  on.exit({
    close(out)
    if (!outComplete) {
      file.remove(destname)
    }
  }, add = TRUE)
  nbytes <- 0L
  repeat {
    bfr <- readBin(inn, what = raw(0L), size = 1L, n = BFR.SIZE)
    n <- length(bfr)
    if (n == 0L) 
      break
    nbytes <- nbytes + n
    writeBin(bfr, con = out, size = 1L)
    bfr <- NULL
  }
  outComplete <- TRUE
  if (remove) {
    close(inn)
    inn <- NULL
    file.remove(filename)
  }
  attr(destname, "nbrOfBytes") <- nbytes
  invisible(destname)
}

## Uncompress and replace file on computer
xzuncompress <- function (filename, destname = gsub("[.]xz$", "", filename, ignore.case = TRUE), 
          temporary = FALSE, skip = FALSE, overwrite = FALSE, remove = TRUE, 
          BFR.SIZE = 1e+07, ...) 
{
  if (!file.exists(filename)) {
    stop("No such file: ", filename)
  }
  if (temporary) {
    destname <- file.path(tempdir(), basename(destname))
  }
  attr(destname, "temporary") <- temporary
  if (filename == destname) {
    stop(sprintf("Argument 'filename' and 'destname' are identical: %s", 
                 filename))
  }
  if (file.exists(destname)) {
    if (skip) {
      return(destname)
    }
    else if (!overwrite) {
      stop(sprintf("File already exists: %s", destname))
    }
  }
  destpath <- dirname(destname)
  if (!file.info(destpath)$isdir) 
    dir.create(destpath)
  inn <- xzfile(filename, open = "rb")
  on.exit(if (!is.null(inn)) close(inn))
  outComplete <- FALSE
  out <- file(destname, open = "wb")
  on.exit({
    close(out)
    if (!outComplete) {
      file.remove(destname)
    }
  }, add = TRUE)
  nbytes <- 0L
  repeat {
    bfr <- readBin(inn, what = raw(0L), size = 1L, n = BFR.SIZE)
    n <- length(bfr)
    if (n == 0L) 
      break
    nbytes <- nbytes + n
    writeBin(bfr, con = out, size = 1L)
    bfr <- NULL
  }
  outComplete <- TRUE
  if (remove) {
    close(inn)
    inn <- NULL
    file.remove(filename)
  }
  attr(destname, "nbrOfBytes") <- nbytes
  invisible(destname)
}
