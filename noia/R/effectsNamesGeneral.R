effectsNamesGeneral <-
function (nloc = 2, max.level = NULL, max.dom = NULL) 
{
    "strrev" <- function(ss) {
        sapply(lapply(strsplit(ss, character(0)), rev), paste, 
            collapse = "")
    }
    ebase <- noia::effectsNames[1:3]
    enames <- ebase
    if (nloc > 1) {
        for (i in 1:(nloc - 1)) {
            enames <- kronecker(ebase, enames, FUN = "paste", 
                sep = "")
            enames <- enames[sapply(enames, statusMaxLevel, max.level)]
            enames <- enames[sapply(enames, statusMaxDom, max.dom)]
        }
    }
    enames <- enames[sapply(enames, statusMaxLevel, max.level)]
    enames <- enames[sapply(enames, statusMaxDom, max.dom)]
    return(strrev(enames))
}
