`varNamesToChar` <-
function (varnam) 
{
    tmp2 <- ""
    tmp1 <- strsplit(varnam, ", ")[[1]]
    for (i in 1:length(tmp1)) {
        tmp2 <- paste(tmp2, tmp1[[i]], "\", \"", sep = "")
    }
    return(tmp2)
}
