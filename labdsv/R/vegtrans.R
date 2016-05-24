vegtrans <- function (taxa,code,value)
{
    if (!is.data.frame(taxa)) {
        taxa <- data.frame(taxa)
    }
    if (length(code) != length(value)) {
        stop("code and value vectors must be of the same length")
    }
    if (is.numeric(code)) {
        code <- c(0,code)
    } else {
        code <- c('0',code)
    }
    if (is.numeric(value)) {
        value <- c(0,value)
    } else { 
        value <- c('0',value)
    }
    newtaxa <- matrix(NA,nrow=nrow(taxa),ncol=ncol(taxa))

    for (i in 1:length(code)) newtaxa[taxa==code[i]] <- value[i]
    newtaxa <- data.frame(newtaxa)
    names(newtaxa) <- names(taxa)
    row.names(newtaxa) <- row.names(taxa)
    if (any(is.na(newtaxa))) {
        print("WARNING, not all values specified")
    }
    return(newtaxa)
}
