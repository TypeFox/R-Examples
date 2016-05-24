"lrzlmomco" <- function(f,para) {
    if(! are.par.valid(para)) return()
    Lu <- vector(mode="numeric", length=length(f))
    MU <- cmlmomco(f=0,para)
    for(i in 1:length(f)) {
       if(f[i] == 0) { Lu[i] <- 0; next }
       tmp <- NULL
       "afunc" <- function(p) {
          return(par2qua(p,para,paracheck=FALSE))
       }
       try(tmp <- integrate(afunc, 0, f[i]))
       Lu[i] <- ifelse(is.null(tmp), NA, tmp$value/MU)
    }
    return(Lu)
}


