mark <-
function(string, view = T, mark = "\""){
    temp <- switch(nchar(mark), rep(mark, 2),
                   c(substring(mark, 1, 1), substring(mark, 2, 2)))
    if(length(temp)) {
        mark <- temp
        temp <- T
    }
    else temp <- F
    if(temp)
        string <- paste(mark[1], paste(string, mark[2], sep = ""), sep = "")
    else warning(paste("Argument \"mark\" should be 1 or 2 characters long",
                       "\n formed by the left and right hand side markers."))
    if(view) cat(string, "\n")
    invisible(string)
}
