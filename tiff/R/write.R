writeTIFF <- function(what, where, bits.per.sample = 8L,
                      compression = c("LZW", "none", "PackBits", "RLE", "JPEG", "deflate"),
                      reduce = TRUE) {
  if (!is.numeric(compression) || length(compression) != 1L) {
    compressions <- c(none=1L, RLE=2L, PackBits=32773L, fax3=3L, fax4=4L, LZW=5L, JPEG=7L, deflate=8L)
    compression <- match.arg(compression)
    compression <- compressions[match(compression, names(compressions))]
  }
  .Call("write_tiff", what, where, bits.per.sample, compression, reduce, PACKAGE="tiff")
}
