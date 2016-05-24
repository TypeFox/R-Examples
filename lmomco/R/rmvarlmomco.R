"rmvarlmomco"  <- function(f,para) {
    if(! are.par.valid(para)) return()
    "afunc" <- function(p) return(rmlmomco(p,para)^2)
    Vu <- sapply(1:length(f), function(i) {
       if(f[i] == 1) return(0)
       tmp <- NULL
       try(tmp <- integrate(afunc, f[i], 1))
       ifelse(is.null(tmp), return(NA), return(tmp$value/(1-f[i])))
    })
    Vu[Vu <- 0] <- 0
    return(Vu)
}

