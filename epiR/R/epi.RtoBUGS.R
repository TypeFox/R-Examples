# Source: Terry Elrod (Terry.Elrod@UAlberta.ca).

"epi.RtoBUGS" <- function(datalist, towhere)
{
    if(!is.list(datalist))
        stop("First argument to writeDatafile must be a list.")
    cat(.formatData(datalist), file = towhere)
}

".formatData" <- function(datalist)
{
    if(!is.list(datalist))
        stop("Argument to formatData must be a list.")
    n <- length(datalist)
    datalist.string <- as.list(rep(NA, n))
    for(i in 1.:n) {
        if(is.numeric(datalist[[i]]) & length(datalist[[i]]) == 1.)
            datalist.string[[i]] <- paste(names(datalist)[i], "=", as.character(datalist[[i]]),
                sep = "")
        if(is.vector(datalist[[i]]) & length(datalist[[i]]) > 1.)
            datalist.string[[i]] <- paste(names(datalist)[i], "=c(", paste(as.character(datalist[[
                i]]), collapse = ","), ")", sep = "")
        if(is.array(datalist[[i]]))
            datalist.string[[i]] <- paste(names(datalist)[i], "=structure(.Data=c(", paste(
                as.character(as.vector(aperm(datalist[[i]]))), collapse = ","), "),.Dim=c(",
                paste(as.character(dim(datalist[[i]])), collapse = ","), "))", sep = "")
    }
    datalist.tofile <- paste("list(", paste(unlist(datalist.string), collapse = ","), ")", sep = "")
    return(datalist.tofile)
}
