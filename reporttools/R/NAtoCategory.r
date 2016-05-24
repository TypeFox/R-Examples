NAtoCategory <- function(fact, label = "missing"){
    if (sum(is.na(fact)) > 0){
        fact <- as.factor(fact)
        leve <- levels(fact)
        leve <- leve[leve != ""]        
        res <- rep(NA, length(fact))
        for (i in 1:length(leve)){res[fact == leve[i]] <- i}
        res[is.na(res) == TRUE] <- length(leve) + 1
        res <- factor(res, levels = 1:(length(leve) + 1), labels = c(leve, 
            label))} else {res <- fact}
    return(res)
}
