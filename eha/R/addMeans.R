addMeans <- function(means, par, var, ns, pfixed, coef.names){
    ## Back-transform the subtracted means;
    ## Will only affect 'scale' and the corresponding row(s)/column(s)
    ## the covariance matrix.

    if (pfixed){
        ncov <- length(par) - ns
        dxy <- diag(ns + ncov)
##        for (i in 1:ns){
        for (i in seq_len(ns)){
            row <- ncov + i
            dxy[row, 1:ncov] <- means
        }
    }else{ # Not pfixed
        ncov <- length(par) - 2 * ns
        dxy <- diag(2 * ns + ncov)
##        for (i in 1:ns){
        for (i in seq_len(ns)){
            row <- ncov + 2 * i - 1
            dxy[row, 1:ncov] <- means
        }
    }

    par <- as.vector(dxy %*% par)
    var <- dxy %*% var %*% t(dxy)

    if (ns > 1){
##        for (i in 1:ns){
        for (i in seq_len(ns)){
            coef.names <- c(coef.names,
                            paste("log(scale)", as.character(i), sep =":"),
                            paste("log(shape)", as.character(i), sep =":"))
        }
        
    }else{
        coef.names <- c(coef.names,
                        "log(scale)", "log(shape)")
    }

    names(par) <- coef.names
    rownames(var) <- colnames(var) <- coef.names
    list(par = par, var = var)
}
