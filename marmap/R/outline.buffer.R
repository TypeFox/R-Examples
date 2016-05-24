outline.buffer <- function(buffer) {
  if (!is(buffer,'bathy')) stop("'buffer' must be a buffer of class bathy")
  buffer[!is.na(buffer)] <- -1
  buffer[is.na(buffer)] <- 0
  return(buffer)
}