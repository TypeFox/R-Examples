"rrmlmomco" <- function(f,para) {
    if(! are.par.valid(para)) return()
    "afunc" <- function(p) return(par2qua(p,para,paracheck=FALSE))
    Qu <- par2qua(f,para,paracheck=FALSE)
    Ru <- sapply(1:length(f), function(i) {
       if(f[i] == 0) return(0)
       tmp <- NULL
       try(tmp <- integrate(afunc, 0, f[i]))
       ifelse(is.null(tmp), return(NA), return(Qu[i] - tmp$value/f[i]))
    })
    Ru[Ru < 0] <- 0 # Truncate to zero, very small negatives might result
    # depending on performance of the numerical integration as well as parameters
    # of the quantile function and the quantile function's performance.
    return(Ru)
}

