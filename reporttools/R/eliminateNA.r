eliminateNA <- function(dat){
    n <- dim(dat)[1]
    c <- dim(dat)[2]
    tmp <- matrix(NA, ncol = c, nrow = n)
    for (i in 1:c){tmp[, i] <- as.numeric(dat[, i])}

    compl <- dat[complete.cases(dat) == TRUE, ]
    incompl <- dat[complete.cases(dat) == FALSE, ]
    res <- list(complete = compl, incomplete = incompl)
    return(res)
}
