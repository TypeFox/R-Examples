###############################################################################
## Read RFLP data
###############################################################################

read.rflp <- function(file){
    x <- read.delim(file = file, stringsAsFactors = FALSE)
    if(!ncol(x) %in% c(3, 4))
        stop("Data in given file", basename(file), "has wrong dimension!")
    
    if(ncol(x) == 3){
        Gel <- sapply(x$Sample, function(x) tail(strsplit(x, "\\_")[[1]], n = 1))
        fun <- function(x){ 
            temp <- strsplit(x, "\\_")[[1]]
            paste(temp[-length(temp)], collapse = "_")
        }
        Sample <- sapply(x$Sample, fun)
        x$Sample <- Sample
        x <- data.frame(x, Gel, stringsAsFactors = FALSE)
    }
    x
}
