#' Brotli Compression
#'
#' Brotli is a compression algorithm optimized for the web, in particular small text
#' documents.
#'
#' Brotli decompression is at least as fast as for gzip while significantly
#' improving the compression ratio. The price we pay is that compression is much
#' slower than gzip. Brotli is therefore most effective for serving static content
#' such as fonts and html pages.
#'
#' For binary (non-text) data, the compression ratio of Brotli usually does not beat
#' \code{bz2} or \code{xz (lzma)}, however decompression for these algorithms is too
#' slow for browsers in e.g. mobile devices.
#'
#' @export
#' @seealso \link{memCompress}
#' @useDynLib brotli R_brotli_compress
#' @name brotli
#' @rdname brotli
#' @param buf raw vector with data to compress/decompress
#' @param mode which compression strategy to use
#' @param quality value between 0 and 11
#' @param log_win log of window size
#' @param log_block log of block size
#' @references J. Alakuijala and Z. Szabadka (October 2015). \emph{Brotli Compressed
#' Data Format}. IETF Internet Draft \url{http://www.ietf.org/id/draft-alakuijala-brotli}.
#' @examples # Simple example
#' myfile <- file.path(R.home(), "COPYING")
#' x <- readBin(myfile, raw(), file.info(myfile)$size)
#' y <- brotli_compress(x)
#' stopifnot(identical(x, brotli_decompress(y)))
#'
#' # Compare to other algorithms
#' length(x)
#' length(brotli_compress(x))
#' length(memCompress(x, "gzip"))
#' length(memCompress(x, "bzip2"))
#' length(memCompress(x, "xz"))
brotli_compress <- function(buf, mode = c("generic", "text", "font"), quality = 11, log_win = 22, log_block = 0){
  stopifnot(is.raw(buf));
  mode <- switch(match.arg(mode),
    generic = 0L,
    text = 1L,
    font = 2L,
    stop("Invalid mode")
  )
  stopifnot(is.numeric(quality))
  stopifnot(is.numeric(log_win))
  stopifnot(is.numeric(log_block))
  .Call(R_brotli_compress, buf, mode, quality, log_win, log_block)
}

#' @export
#' @useDynLib brotli R_brotli_decompress
#' @rdname brotli
brotli_decompress <- function(buf){
  stopifnot(is.raw(buf))
  .Call(R_brotli_decompress, buf)
}
