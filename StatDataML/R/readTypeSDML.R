readType <- function(x)
{
    ret <- list()
    ret$type <- "character"
    ret$mode <- NULL
    ret$labels <- NULL
    ret$min <- -Inf
    ret$max <- +Inf

    if (xmlName(x) == "type") {
        ## FIXME: why is the seq_along() statement needed??
        x <- x[sapply(x[seq_along(x)],
                      xmlName) != "text"][[1]]
        ret$type <- xmlName(x)

        if (ret$type == "categorical") {
            ## mode
            ret$mode <- xmlAttrs(x)["mode"]

            ## labels
            ind <- 1
            for (k in 1:xmlSize(x)) {
                if (xmlName(x[[k]]) != "text") {
                    code <- as.integer(xmlAttrs(x[[k]])["code"])
                    children <- xmlChildren(x[[k]])
                    lab <- if (length(children)) children else NULL
                    ret$labels[code] <-
                        if (is.null(lab)) ind else xmlValue(lab[[1]])
                    ind <- ind + 1
                }
            }
        }

        if (ret$type == "numeric") {
            if (!length(xmlChildren(x)))
                ret$mode <- "real"
            else {
                x <- x[sapply(x[seq_along(x)],
                              xmlName) != "text"][[1]]
                ret$mode <- xmlName(x)
                l <- length(xmlChildren(x))
                if (ret$mode %in% c("integer", "real") && l)
                    for (i in seq_len(l))
                        if (xmlName(x[[i]]) != "text")
                            ret[[xmlName(x[[i]])]] <- getDataSDML(x[i])
            }
        }
    }

    ret
}
