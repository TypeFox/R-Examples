#  -----------------------------------------------------------
#  pySource
#  ========
## @title Reads mixed R and Python code from a file
##
## @description The function BEGIN.Python allows interactive development
##              but doesn't work in combination with the function source.
##              Therefore pySource provides an alternative to the function source
##              which also can handle BEGIN.Python statements.
## @param file a character string giving the pathname of the file.
## @param local see documentation of source
## @param echo see documentation of source
## @param print.eval see documentation of source
## @param verbose see documentation of source
## @param prompt.echo see documentation of source
## @param max.deparse.length see documentation of source
## @param chdir see documentation of source
## @param encoding see documentation of source
## @param continue.echo see documentation of source
## @param skip.echo see documentation of source
## @param keep.source see documentation of source
## @details The function pySource works similar to source, but code 
##          which is enclosed between BEGIN.Python and END.Python is
##          replaced by pyExec and the quoted version of the code.
## @examples
## \dontrun{
## writeLines(c("x <- 3", "BEGIN.Python()", 
##              "x=3**3", "print(3*u'Hello R!\\n')", 
##              "END.Python"), "myMixedCode.R")
## pySource("myMixedCode.R")
## }
#  -----------------------------------------------------------
pySource <- function(file, local = FALSE, echo = verbose, print.eval = echo, 
    verbose = getOption("verbose"), prompt.echo = getOption("prompt"), 
    max.deparse.length = 150, chdir = FALSE, encoding = getOption("encoding"), 
    continue.echo = getOption("continue"), skip.echo = 0, 
    keep.source = getOption("keep.source")){
    code <- readLines(file)
    temp <- tempfile()
    code <- paste(code, collapse="\n")

    m <- unlist(regmatches(code,
                           gregexpr("BEGIN\\.Python\\(\\)(.*?)END\\.Python", code)))
    repl <- gsub("END.Python", "", gsub("BEGIN.Python()", "", m, fixed=T), fixed=T)
    repl <- sprintf("pyExec(%s)", shQuote(repl))
    for (i in 1:length(m)){
        code <- sub(m[i], repl[i], code, fixed=TRUE)
    }

    writeLines(code, temp)
    source(temp, local, echo, print.eval, verbose, prompt.echo, 
           max.deparse.length, chdir, encoding, continue.echo, skip.echo,
           keep.source)
}

#  -----------------------------------------------------------
#  BEGIN.Python
#  ============
#' @title Execute Python interactively from within R
#'
#' @description The function BEGIN.Python starts an Python read-eval-print loop.
#' @details BEGIN.Python emulates the behavior of the Python terminal
#'          and therefore allows interactive Python code development
#'          from within R.
#' @return Returns the entered code as character, code lines which throw an
#'         exception are omitted.
#' @note This won't work with RStudio because of a known  
#' \href{https://support.rstudio.com/hc/communities/public/questions/206744317-readLines-Bug-?page=1#answer-206647788}{RStudio issue}.
#' @examples
#' \dontshow{PythonInR:::pyCranConnect()}
#' \dontrun{
#' code <-
#' BEGIN.Python()
#' import os
#' os.getcwd()
#' dir(os)
#' x = 3**3
#' for i in xrange(10):
#'     if (i > 5):
#'         print(i)
#'
#' END.Python
#' ## NOTE: BEGIN.Python returns the successfully executed code as character.
#' cat(code, sep="\n")
#' pyGet0("x")
#' }
#  -----------------------------------------------------------
BEGIN.Python <- function(){
    if ( pyConnectionCheck() ) return(invisible(NULL))
    f <- file("stdin")
    open(f)
    cat("py> ")
    pyCode <- character()
    execBuffer <- ""
    while(TRUE) {
        line <- readLines(f,n=1)
        if (grepl("END.Python", line, fixed=TRUE)) {
            if ( nchar(execBuffer) > 0 ) pyExec(execBuffer)
            break
        }
        if (grepl("(^\\s|:\\s*$)", line)) {
            execBuffer <- paste(c(execBuffer, line), collapse="\n")
        } else {
            if ( nchar(execBuffer) > 0 ) {
                pyExec(execBuffer)
                cat("py> ")
            }
            if ( nchar(line) > 0 ) {
            	tryCatch({
            		pyExecp(line)
                	pyCode <- c(pyCode, line)
                  	},
                 	warning=function(w){ print(w) },
                 	error=function(e){ print(e) }
                 	)
            	cat("py> ")
        	}
        }
    }
    return(invisible(pyCode))
}

