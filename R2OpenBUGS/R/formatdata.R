"formatdata" <-
function (datalist){
    if (!is.list(datalist) || is.data.frame(datalist)) 
        stop("Argument to formatdata() must be a list.")
    n <- length(datalist)
    datalist.string <- vector(n, mode = "list")
    for (i in 1:n) {
        if (length(datalist[[i]]) == 1) 
            datalist.string[[i]] <- paste(names(datalist)[i], 
                "=", as.character(datalist[[i]]), sep = "")
        if (is.vector(datalist[[i]]) && length(datalist[[i]]) > 1) 
            datalist.string[[i]] <- paste(names(datalist)[i], 
                "=c(", paste(as.character(datalist[[i]]), collapse = ", "), 
                ")", sep = "")
        if (is.array(datalist[[i]])) 
            datalist.string[[i]] <- paste(names(datalist)[i], 
                "= structure(.Data= c(", paste(as.character(as.vector(aperm(datalist[[i]]))), 
                  collapse = ", "), "), .Dim=c(", paste(as.character(dim(datalist[[i]])), 
                  collapse = ", "), "))", sep = "")
    }
    datalist.tofile <- paste("list(", paste(unlist(datalist.string), 
        collapse = ", "), ")", sep = "")
    datalist.tofile
}
