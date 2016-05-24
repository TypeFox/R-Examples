Capitalize <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
}

FormatURL <- function(paths, query) {
    pathstr <- paste(paths, collapse="/")
    querystr <- paste(names(query), query, sep="=", collapse="&")
    url <- paste(pathstr, querystr, sep="?")
    return (url)
}

IsValidReply <- function(reply) {
    return(is.element(reply, c("y", "n")))
}
