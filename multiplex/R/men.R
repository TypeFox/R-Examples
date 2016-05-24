men <-
function (y, prsep = ", ") 
{
    mn <- vector()
    if (isTRUE(length(y) > 0) == TRUE) 
        for (i in 1:length(y)) if (isTRUE(strsplit(y[i], prsep)[[1]][1] < 
            strsplit(y[i], prsep)[[1]][2]) == TRUE) 
            mn <- append(mn, y[i])
    rm(i)
    return(mn)
}
