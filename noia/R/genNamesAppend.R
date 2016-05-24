genNamesAppend <-
function (name) 
{
    ans <- NULL
    for (n in name) {
        for (g in noia::genotypesNames[1:3]) {
            ans <- c(ans, paste(g, n, sep = ""))
        }
    }
    return(ans)
}
