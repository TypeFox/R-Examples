getData <- function(path)
{
    files <- c("anneal.data", "anneal.test")
    files <- paste(path, files, sep="/")
    dat1 <- read.table(files[1],
        sep=",",
        comment.char="",
        na.strings="",
        stringsAsFactors=FALSE, strip.white=TRUE)
    dat2 <- read.table(files[2],
        sep=",",
        comment.char="",
        na.strings="",
        stringsAsFactors=FALSE, strip.white=TRUE)
    dat <- rbind(dat1, dat2)
    rm(dat1, dat2)
    for (i in 1:ncol(dat)) {
        if (class(dat[, i]) == "character") {
            dat[dat[, i] == "?", i] <- "-"
        }
    }
    dat
}

