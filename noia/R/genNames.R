genNames <-
function (nloc = 2) 
{
    names <- noia::genotypesNames[1:3]
    if (nloc > 1) {
        for (i in 2:nloc) {
            names <- genNamesAppend(names)
        }
    }
    return(names)
}
