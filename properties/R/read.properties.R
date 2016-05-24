#' Read a Java Properties File for Use in RSB Applications
#' @param file properties file, either a character string (path) or a connection 
#' @param fields subset of field names, if NULL, all fields are included
#' @param encoding encoding of the properties file to read (default value \code{"UTF-8"})
#' @return list with key-value pairs 
#' @author Tobias Verbeke
#' @export
read.properties <- function(file, fields = NULL, encoding = "UTF-8"){
    
    if (is.character(file)) {
        file <- file(file, "r", encoding = encoding)
        on.exit(close(file))
    }
    if (!inherits(file, "connection")) 
        stop("'file' must be a character string or connection")
    
    lines <- readLines(file)
    commentedLines <- grepl("^#.*$", lines)
    lines <- lines[!commentedLines]
    
    line_is_not_empty <- !grepl("^[[:space:]]*$", lines)
    lines <- lines[line_is_not_empty]

    line_has_tag <- grepl("^[^[:blank:]][^=:]*=", lines)
    ind <- which(!line_has_tag)
    
    if (length(ind)) {
        lines <- strtrim(lines[ind], 0.7 * getOption("width"))
        stop("Invalid Java properties file format.\nContinuation lines must not start a record.\nOffending lines start with:\n%s", 
             paste("  ", lines, sep = "", collapse = "\n"))
    }
    keys <- gsub("^([^=]+)=.*$", "\\1", lines)
    values <- gsub("^[^=]+=(.*)$", "\\1", lines)
    
    names(values) <- keys
    out <- as.list(values)
    out <- if (!is.null(fields)) out[fields] else out
    return(out)
} 
