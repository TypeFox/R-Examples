
HSAURcite <- function(pkg) {
    ct <- citation(pkg)
    attr(ct, "label") <- paste("PKG:", pkg, sep = "", collapse = "")
    for (n in c("note"))
        ct[[n]] <- gsub("R", "\\R{}", ct[[n]])
    class(ct) <- "HSAURcitation"
    return(ct)
}

toBibtex.HSAURcitation <-  function (object, ...) 
{
    z <- paste("@", attr(object, "entry"), "{", attr(object, "label"), 
               ",", sep = "")
    if ("author" %in% names(object)) {
        object$author <- toBibtex(object$author)
    }
    for (n in names(object)) z <- c(z, paste("  ", n, " = {", 
        object[[n]], "},", sep = ""))
    z <- c(z, "}")
    class(z) <- "Bibtex"
    z
}
