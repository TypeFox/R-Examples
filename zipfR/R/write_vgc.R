write.vgc <- function (vgc, file)
{
  if (! inherits(vgc, "vgc")) stop("'vgc' argument must be of class 'vgc'")
  if ( (length(file) != 1) || (! is.character(file)) )
    stop("'file' argument must be a single character string")
  
  write.table(vgc, file=auto.gzfile(file), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}
