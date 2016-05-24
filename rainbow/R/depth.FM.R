`depth.FM` <- function(data, trim = 0.25, xeps = 0.00000001, x = NULL){
    functions = t(data$y)
    nrow <- dim(functions)[1]
    ncol <- dim(functions)[2]
    if(is.null(nrow) && is.null(ncol))
       stop("I do not have a matrix")
    if(is.null(x)) 
       x = 1:ncol
    d <- matrix(NA, nrow = nrow, ncol = ncol)
    for(i in 1:ncol){
        Fn = ecdf(functions[,i])
        d[,i] = 1 - abs(0.5 - Fn(functions[,i]))
    }
    ans <- apply(d, 1, mean)
    k = which.max(ans)
    med = functions[k,]
    lista = which(ans >= quantile(ans, probs = trim, na.rm = TRUE))
    mtrim = apply(functions[lista,], 2, mean)
    return(list("median" = med, "lmed" = k, "mtrim" = mtrim, 
           "ltrim" = lista, "prof" = ans))
}

