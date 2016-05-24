"rrmvarlmomco"  <- function(f,para) {
    if(! are.par.valid(para)) return()
    "afunc" <- function(p) return(rrmlmomco(p,para)^2)
    Du <- sapply(1:length(f), function(i) {
       if(f[i] == 0) return(0)
       tmp <- NULL
       try(tmp <- integrate(afunc, 0, f[i]))
       ifelse(is.null(tmp), return(NA), return(tmp$value/f[i]))
    })
    Du[Du < 0] <- 0
    return(Du)
}

