"pdfgep" <-
function(x, para) {
    if(! are.pargep.valid(para)) return()
    names(para$para) <- NULL
    B <- 1/para$para[1]
    K <-   para$para[2]
    H <-   para$para[3]

    x[x < 0] <- NA
    HEX <- H*exp(-B*x)
    e1 <- exp(-H +       HEX)
    e2 <- exp(-H - B*x + HEX)
    C  <- (K*B*H)/(1-exp(-H))^K
    f <- C*(1 - e1)^(K - 1)*e2

    names(f) <- NULL
    f[! is.finite(f)] <- NA
    f[is.na(f)] <- 0 # decision Dec. 2015
    return(f)
}

