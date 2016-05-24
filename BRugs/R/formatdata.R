"formatdata" <-
function (datalist){
    if (!is.list(datalist) || is.data.frame(datalist))
      stop("argument to formatdata() ", "must be a list")
    n <- length(datalist)
    datalist.string <- vector(n, mode = "list")
    datanames <- names(datalist)
    for (i in 1:n) {
        if (is.factor(datalist[[i]]))
            datalist[[i]] <- as.integer(datalist[[i]])
        datalist.string[[i]] <-
        if (length(datalist[[i]]) == 1)
            paste(names(datalist)[i],
                  "=", as.character(datalist[[i]]), sep = "")
        else if (is.vector(datalist[[i]]) && length(datalist[[i]]) > 1)
            paste(names(datalist)[i],
                  "=c(", paste(as.character(datalist[[i]]), collapse = ", "),
                  ")", sep = "")
        else
            paste(names(datalist)[i],
                  "= structure(.Data= c(",
                  paste(as.character(as.vector(aperm(datalist[[i]]))), collapse = ", "),
                  "), .Dim=c(",
                  paste(as.character(dim(datalist[[i]])), collapse = ", "),
                  "))", sep = "")
    }
    datalist.tofile <- paste("list(",
        paste(unlist(datalist.string), collapse = ", "),
        ")", sep = "")
    return(datalist.tofile)
}
