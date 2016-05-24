"quakap" <-
function(f,para,paracheck=TRUE) {
    if(! check.fs(f)) return()
    if(paracheck == TRUE) {
      if(! are.parkap.valid(para)) return()
    }
    U <- para$para[1]
    A <- para$para[2]
    G <- para$para[3]
    H <- para$para[4]

    x <- sapply(1:length(f), function(i) {
                    if(f[i] == 0) {
                       if(H <= 0 & G  < 0) return(U+A/G)
                       if(H >  0 & G != 0) return(U+A/G*(1-H^-G))
                       if(H >  0 & G == 0) return(U+A*log(H))
                       if(H <= 0 & G >= 0) return(-Inf)
                       stop("f is fine: should not be here in code execution.")
                    } else if(f[i] == 1) {
                       if(G <= 0) return(Inf)
                       return(U+A/G)
                       stop("f=1: should not be here in code execution.")
                    } else {
                       Y <- -log(f[i])
                       if(H != 0) Y <- (1-exp(-H*Y))/H
                       Y <- -log(Y)
                       if(G != 0) Y <- (1-exp(-G*Y))/G
                       return(U+A*Y)
                    } })
    names(x) <- NULL
    return(x)
}

