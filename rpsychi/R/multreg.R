multreg <- function(formula, data, sig.level=.05, digits=3){
    ##setting
    depname <- unlist(strsplit(as.character(formula), " ")[2])
    x <- unlist(strsplit(as.character(formula), " ")[3])
    indname <- x[!is.element(x, "+")]
    x <- c(depname, indname)
           
    ##sample statistics
    samp.stat <- data.frame(m = apply(data[, x], 2, mean), sd = apply(data[, x], 2, sd))
    
    r <- cor(data[, x])[1,2]
    n <- nrow(data)
    
    

    output <- multreg.second(corr=cor(data[,x]), formula=formula, n=nrow(data), 
    m = apply(data[, x], 2, mean), sd = apply(data[, x], 2, sd), sig.level=sig.level, digits=digits)
    output <- sapply(c(list(samp.stat=samp.stat), output), round, digits)
    return(output)
}
