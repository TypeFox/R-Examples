"rmlmomco" <- function(f,para) {
    if(! are.par.valid(para)) return()
    "afunc" <- function(p) return(par2qua(p,para,paracheck=FALSE))
    Qu <- par2qua(f,para,paracheck=FALSE)
    Mu <- sapply(1:length(f), function(i) {
       if(f[i] == 1) return(0)
       tmp <- NULL
       try(tmp <- integrate(afunc, f[i], 1))
       ifelse(is.null(tmp), return(NA), return(tmp$value/(1-f[i]) - Qu[i]))
    })
    Mu[Mu < 0] <- 0 # At least for the example in docs, as the residual
    # life goes to zero via the integration and deep into the rigth tail,
    # numerical variations touching just below zero, say -1E-12 can happen,
    # so truncate the result.
    return(Mu)
}
