as.data.frame.registry <-
function(x, ...)
{
    .one_line <- function(entry) {
        entry <- lapply(.labels(entry), function(i) i[[1]])
        data.frame(unclass(entry), ...)
    }

    ret <- do.call(rbind, lapply(x$get_entries(), .one_line))
    row.names(ret) <- NULL
    ret
}

.labels <-
function(x)
{
    ## transform function entries
    x[sapply(x, inherits, "function")] <- "<<function>>"

    ## transform objects
    obj <- sapply(x, is.object)
    x[obj] <- paste("<<", sapply(x, class)[obj], ">>", sep = "")

    x
}
