## a constructor for named lists
nlist <- function(...) {
    x <- list(...)
    names(x) <- character(0)
    x
}

## ----------------------------------------------------------------------------- 
##
##   Not Exported Function for Testing
##
## -----------------------------------------------------------------------------
sys_call <- function(cmd, args) {
    x <- tryCatch(system2(cmd, args, stdout=TRUE, stderr=TRUE), error=function(e) NULL)
    return(x)
}

pandoc_to_json <- function(file, from="markdown") {
    args <- sprintf("-f %s -t json %s", from, file)
    system2("pandoc", args, stdout=TRUE, stderr=TRUE)
}

pandoc_from_json <- function(json, to) {
    args <- sprintf("%s | pandoc -f json -t %s", shQuote(json), to)
    system2("echo", args, stdout=TRUE, stderr=TRUE)
}

#  -----------------------------------------------------------
#  pandocfilters_writer
#  ====================
#' @title Write the JSON-formatted AST to a connection
#' @description Write the JSON-formatted AST to a connection.
#' @param x a JSON representation of the AST to be written out
#' @param con a connection object or a character string to which the JSON-formatted AST is written
#' @param format a character string giving the format (e.g. \code{"latex"}, \code{"html"})
#' @details If you want to apply a filter to the document before it get's written out, or your
#'          pandoc installation is not registered in the \code{PATH} it can be favorable to provide your
#'          own writer function to the document class.
#' @export
#  -----------------------------------------------------------
pandocfilters_writer <- function(x, con, format) {
    args <- sprintf("%s | pandoc -f json -t %s", shQuote(as.character(x)), format)
    x <- system2("echo", args, stdout=TRUE, stderr=TRUE)
    writeLines(x, con=con)
}

test <- function(x, to="html") {
    d <- list(list(unMeta=nlist()), x)
    pandoc_from_json(as.character(jsonlite::toJSON(d, auto_unbox=TRUE)), to=to)
}

detect_pandoc_version <- function() {
    x <- sys_call("pandoc", "--version")
    if ( is.null(x) ) {
        writeLines("\n\nInfo message:\nCouldn't find 'pandoc'!\nPandoc version is set to '1.16', use `set_pandoc_version` to change the settings if necessary.\n\n")
        return(NULL)
    }
    b <- grepl("pandoc +\\d\\.\\d", x)  
    if ( !any(b) ) return(NULL)
    x <- x[b][1]
    version <- gsub("[[:alpha:] ]", "", x)
    version_numeric <- as.numeric(regmatches(version, regexpr("\\d+\\.\\d+", version)))
    list(char=version, num=version_numeric)
}
## detect_pandoc_version()

#  -----------------------------------------------------------
#  get_pandoc_version
#  ==================
#' @title Get Pandoc Version
#' @description Get the version of pandoc.
#' @export
#  -----------------------------------------------------------
get_pandoc_version <- function() {
    base::getNamespace("pandocfilters")$pandoc$version
}

#  -----------------------------------------------------------
#  set_pandoc_version
#  ==================
#' @title Set Pandoc Version
#' @description Set the version version pandoc.
#' @param x a numeric giving the pandoc version (e.g. 1.14 or 1.15 or 1.16 or 1.17)
#' @export
#  -----------------------------------------------------------
set_pandoc_version <- function(x) {    
    assign("version", x, envir=base::getNamespace("pandocfilters")$pandoc)
}
