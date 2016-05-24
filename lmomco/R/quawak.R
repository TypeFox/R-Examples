"quawak" <-
function(f,wakpara,paracheck=TRUE) {
    if(! check.fs(f)) return()
    if(paracheck == TRUE) {
      if(! are.parwak.valid(wakpara)) return()
    }
    #    #   UFL SHOULD BE CHOSEN SO THAT EXP(UFL) JUST DOES NOT CAUSE
    #    UNDERFLOW
    #
    UFL <- log(.Machine$double.xmin);
    XI <- wakpara$para[1]
    A <- wakpara$para[2] # alpha
    B <- wakpara$para[3] # beta
    C <- wakpara$para[4] # gamma
    D <- wakpara$para[5] # delta

    x <- sapply(1:length(f), function(i) {
                       if(f[i] == 0) return(XI)
                       if(f[i] == 1) return(XI+A/B*(1-0^B) - C/D*(1 - 0^(-D)))
                       Y1 <- Z <- -log(1-f[i])
                       if(B == 0) {
                          Y2 <- Z
                          if(D != 0) Y2 <- (1-exp(D*Y2))/(-D)
                          return(XI+A*Y1+C*Y2)
                       } else {
                          TEMP <- -B*Z
                          if(TEMP <  UFL) Y1 <- 1/B
                          if(TEMP >= UFL) Y1 <- (1-exp(TEMP))/B
                          Y2 <- Z
                          if(D != 0) Y2 <- (1-exp(D*Y2))/(-D)
                          return(XI+A*Y1+C*Y2)
                       }
                     })
    names(x) <- NULL
    return(x)
}

