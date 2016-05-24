zero.r.test <- function(formula, data, sig.level=.05, digits=3){

    ##setting
    depname <- unlist(strsplit(as.character(formula), " ")[2])
    indname <- unlist(strsplit(as.character(formula), " ")[3])
    x <- c(depname, indname)
    
    ##sample statistics
    samp.stat <- data.frame(m = apply(data[, x], 2, mean), sd = apply(data[, x], 2, sd))
    
    r <- cor(data[, x])[1,2]
    n <- nrow(data)

    output <- zero.r.test.second(r=r, n=n, sig.level=sig.level, digits=digits)
    output <- sapply(c(list(samp.stat=samp.stat), output), round, digits)
    return(output)
}
