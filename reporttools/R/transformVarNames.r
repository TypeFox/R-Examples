transformVarNames <- function(dat, name){
    nams <- names(dat)
    if (is.null(nams) == TRUE){nams <- dimnames(dat)[[1]]}
    ns <- nchar(nams)
    out <- rep(NA, length(ns))
    for (i in 1:length(ns)){
        space <- ""
        for (j in 1:(max(ns) - ns[i] + 1)){
            space <- paste(space, " ", sep = "")
        }
        out[i] <- paste(nams[i], space, '<- ', name, '[ ,"', nams[i], '"]', sep = "")
    }
    return(cat(matrix(out, ncol = 1), sep = "\n"))
}
