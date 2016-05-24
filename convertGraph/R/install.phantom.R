#' @title Saves Absolute Path to PhantomJS
#' @keywords PhantomJS
#' @usage install.phantom(path)
#' @return The function does not return anything
#' @description takes the path to executable
#' \code{'phantomJS'} binary and memorizes its absolute path for future calls. PhantomJS binary
#' can be downloaded from \url{http://phantomjs.org/download.html}
#' @param path Path to executable 'phantomJS' binary, downloadble from \url{http://phantomjs.org/download.html}
#'
#' @author E. F. Haghish \cr
#' Medical Informatics and Biostatistics (IMBI) \cr
#' University of Freiburg, Germany \cr
#' \email{haghish@imbi.uni-freiburg.de} \cr
#' \cr
#' Department of Mathematics and Computer Science \cr
#' University of Southern Denmark \cr
#' \email{haghish@imada.sdu.dk}
#'
#'
#' @examples
#' \dontrun{
#' #save the absolute path to phantomJS by giving relative path
#' install.phantom("./bin/phantomjs")
#'
#' #save the absolute path to phantomJS
#' install.phantom("c:\phantomjs\bin\phantomjs")
#'
#' #convert JPEG to PDF
#' convertGraph("./example.jpeg", "./example.pdf", path = "path to executable phantomJS" )
#' }
#'
#' @export
#' @importFrom tools file_path_as_absolute
#' @importFrom tools file_ext

install.phantom <- function(path) {
    if (file.exists(path)) {
        path <- file_path_as_absolute(path)
        memory <- system.file("PATHMEMORY", package = "convertGraph")
        write(path, file=memory, append=TRUE)
    }
}
