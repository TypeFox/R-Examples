
file.head <- function(file, n = 6, truncate.cols = TRUE){
    lns <- readLines(file, n = n)
    lns <- gsub("\t", "\\\\t", lns)
    if(truncate.cols){
        try(file.head.lns <- substr(lns, 1, getOption("width") - 1), silent = TRUE)
        if(!exists("file.head.lns")){
            Encoding(lns) <- "bytes"
            file.head.lns <- substr(lns, 1, getOption("width") - 1)
        }
    }
    cat(file.head.lns, sep = "\n")
}

