mmap.csv <-
function (file, header = TRUE, sep = ",", quote = "\"", dec = ".", 
    fill = TRUE, comment.char = "", row.names, ...) 
{
  ncols <- length(gregexpr(sep,readLines(file,1))[[1]]) + 1
  mcsv <- tempfile()
  tmplist <- vector("list",ncols)
  cnames <- character(ncols)
  if( !missing(row.names) && is.numeric(row.names) && length(row.names)==1L)
    ncols <- ncols-1
  for(col in 1:ncols) {
    colclasses <- rep("NULL",ncols)
    colclasses[col] <- NA
    clm <- read.table(file=file, header=header, sep=sep, quote=quote, dec=dec,
                      fill=fill, comment.char=comment.char, colClasses=colclasses, stringsAsFactors=FALSE,
                      row.names=row.names,...)
    cnames[col] <- colnames(clm)
    tmplist[[col]] <- as.mmap(clm[,1], force=TRUE) 
  }
  stype <- do.call(struct,lapply(tmplist, function(X) X$storage.mode))
  totalsize <- sum(sapply(tmplist, nbytes))
  tmpstruct <- tempfile()
  writeBin(raw(totalsize), tmpstruct)
  tmpstruct <- mmap(tmpstruct, stype)
  for(col in 1:ncols) {
    tmpstruct[,col] <- tmplist[[col]][]
  }
  colnames(tmpstruct) <- cnames
  extractFUN(tmpstruct) <- as.data.frame
  tmpstruct
}

