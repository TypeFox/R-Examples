swp <-
function (y) 
{
    sy <- vector()
    length(sy) <- length(y)
    for (i in 1:length(y)) sy[i] <- paste(strsplit(y[i], ", ")[[1]][2], 
        strsplit(y[i], ", ")[[1]][1], sep = ", ")
    rm(i)
    return(sy)
}
