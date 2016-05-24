"quape3" <-
function(f,para,paracheck=TRUE) {
    if(! check.fs(f)) return()
    if(paracheck == TRUE) {
      if(! are.parpe3.valid(para)) return()
    }
    # SMALL IS USED TO TEST WHETHER SKEWNESS IS EFFECTIVELY ZERO
    SMALL <- sqrt(.Machine$double.eps)
    U <- para$para[1]
    A <- para$para[2]
    G <- para$para[3]

    x <- vector(mode="numeric", length=length(f))
    if(abs(G) <= SMALL) {
       x <- U+A*qnorm(f)
    } else {
       ALPHA <- 4/G^2
       BETA <- abs(0.5*A*G)
       if(G > 0) {
          x <- U-ALPHA*BETA+qgamma(  f,ALPHA,scale=BETA)
       } else {
          x <- U+ALPHA*BETA-qgamma(1-f,ALPHA,scale=BETA)
       }
    }

    x[f == 0 & G > 0] <- U-2*A/G
    x[f == 1 & G < 0] <- U-2*A/G

    names(x) <- NULL
    return(x)
}
