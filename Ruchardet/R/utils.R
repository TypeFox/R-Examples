#' detect file encoding for inputted file 
#' @param file to read from, or a \code{\link{connection}} which will be opened if necessary, and if so closed at the end of the function call.
#' @param n integer. The (maximal) number of lines to read. Negative values indicate that one should read up to the end of input on the connection.
#' @param default default encoding if fail to expect.
#' @return encoding name
#' @examples 
#' big5encfile <- file.path(system.file(package="Ruchardet"),"tests","big5.txt") 
#' detectFileEncoding(big5encfile) 
#' \dontrun{
#' detectFileEncoding("http://www.ppomppu.co.kr/")
#' detectFileEncoding("http://freesearch.pe.kr")
#' }
#' @export 
detectFileEncoding <- function(file, n=-1, default=getOption("encoding")){

  if (is.character(file)) {
      file <- file(file, "rt")
      on.exit(close(file))
  }
  if (!inherits(file, "connection")) 
      stop("'file' must be a character string or connection")
  if (!isOpen(file, "rt")) {
      open(file, "rt")
      on.exit(close(file))
  }
  encs <- detectEncoding(paste0(readLines(file, warn=FALSE, n=n),collapse=""))
  if(encs == ''){
    warning("can't expect encoding, will return 'default' encoding")
    encs <- default
  }
  return(encs)
}




