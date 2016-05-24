##' Get CXXFLAGS to compile with covafill
##'
##' @title CXXFLAGS to compile with covafill
##' @return Returns a string with the CXXFLAGS needed to compile C++ code using covafill.
##' @author Christoffer Moesgaard Albertsen
##' @seealso \code{\link[TMB]{compile}}
##' @examples
##' \dontrun{
##' if(require("TMB")){
##'    f <- system.file("examples","tmbtest","tmbtest.cpp", package='covafillr')
##'    TMB::compile(f,CXXFLAGS = cxxFlags())
##' }
##' }
##' @export
cxxFlags <- function(){
    paste0("-I",system.file("include",package="covafillr"))
}
