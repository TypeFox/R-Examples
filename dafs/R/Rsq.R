Rsq = function(R, fmt = "$R^2 = %6.4f$", adj = FALSE){
    if(adj){
        fmt = paste("adjusted", fmt)
    }
    return(sprintf(fmt, R))
}
