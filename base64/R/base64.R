#' Encode and Decode base64
#'
#' Wraps \code{\link[openssl:base64_encode]{openssl::base64_encode}}
#' to replace the deprecated implementation by Romain Francois.
#'
#' @export
#' @rdname base64
#' @aliases base64
#' @param input input file
#' @param output output file
#' @param linebreaks insert linebreaks to make output human readable
#' See \code{\link[openssl:base64_encode]{openssl::base64_encode}}
#' @importFrom openssl base64_encode
#' @examples # encode a file
#' myfile <- R.home("COPYING")
#' tmp <- tempfile()
#' base64::encode(myfile, tmp)
#'
#' # decode it back
#' orig <- tempfile()
#' base64::decode(tmp, orig)
#' readLines(orig)
encode <- function(input, output = tempfile(), linebreaks = TRUE){
  input <- normalizePath(input, mustWork = TRUE)
  buf <- readBin(input, raw(), file.info(input)$size)
  base64 <- base64_encode(buf, linebreaks = linebreaks)
  writeLines(base64, output)
  output
}

#' @export
#' @rdname base64
#' @importFrom openssl base64_decode
decode <- function(input, output = tempfile()){
  input <- normalizePath(input, mustWork = TRUE)
  buf <- readBin(input, raw(), file.info(input)$size)
  bin <- base64_decode(buf)
  writeBin(bin, output)
  output
}

#' Encode a png file as a img data uri
#'
#' This creates html code to embed a png file into an html document.
#' \Sexpr[results=rd, stage=build, echo=FALSE]{
#'   pngfile <- tempfile()
#'   png(pngfile, width = 600, height = 400 )
#'   plot(1:100, rnorm(100), pch = 21, bg = "red", cex = 2 )
#'   dev.off()
#'   base64::img( pngfile, Rd = TRUE )
#'}
#'
#' @export
#' @importFrom openssl base64_encode
#' @param file png file to translate into a data uri
#' @param Rd if \code{TRUE}, extra markup is added to facilitate inclusion
#' of the image in an Rd file
#' @param alt alternate text
#' @examples pngfile <- tempfile()
#' png(pngfile, width = 600, height = 400)
#' plot(1:100, rnorm(100), pch = 21, bg = "red", cex = 2 )
#' dev.off()
#' img(pngfile, Rd = TRUE)
img <- function( file, Rd = FALSE, alt = "image" ){
  input <- normalizePath(file, mustWork = TRUE)
  buf <- readBin(input, raw(), file.info(input)$size)
  base64 <- base64_encode(buf, linebreaks = FALSE)
  sprintf('%s<img src="data:image/png;base64,\n%s" alt="%s" />%s',
    if( Rd ) "\\out{" else "", base64, alt, if( Rd ) "}" else ""
  )
}
