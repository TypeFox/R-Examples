write.spc <- function (spc, file)
{
  if (! inherits(spc, "spc")) stop("'spc' argument must be of class 'spc'")
  if ( (length(file) != 1) || (! is.character(file)) )
    stop("'file' argument must be a single character string")
  if (attr(spc, "m.max") > 0)
    warning("saving incomplete frequency spectrum, which cannot be restored from disk file!")
  if (attr(spc, "hasVariances"))
    warning("variance of expected vocabulary size cannot be saved to disk file!")
  
  write.table(spc, file=auto.gzfile(file), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}
