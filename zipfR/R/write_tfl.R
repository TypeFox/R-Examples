write.tfl <- function (tfl, file)
{
  if (! inherits(tfl, "tfl")) stop("'tfl' argument must be of class 'tfl'")
  if ( (length(file) != 1) || (! is.character(file)) )
    stop("'file' argument must be a single character string")
  if (attr(tfl, "incomplete"))
    warning("saving incomplete type frequency list, which cannot be restored from disk file!")
  
  write.table(tfl, file=auto.gzfile(file), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}
