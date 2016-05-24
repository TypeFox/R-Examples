sdfhd <- function(n, alpha=1, phi=numeric(0), theta=numeric(0), lmodel=c("FD", "FGN", "PLA", "NONE")) {
    lmodel <- match.arg(lmodel)
    if(!(InvertibleQ(phi)&&InvertibleQ(theta)&&alpha>0&&alpha<2)) 
        stop("error: non-invertible or non-stationary")
    if(length(phi)==0&&length(theta)==0) s1 <- 1 else s1 <- 2*pi*sdfarma(n, phi, theta)
    s2 <- switch(lmodel,
     FD   = sdfFD((1-alpha)/2, n),
     FGN  = sdfFGN(1-alpha/2, n),
     PLA  = sdfPLA(alpha, n),
     PLS  = sdfPLS(alpha, n),
     NONE = sdfFGN(0.5, n))
     s1*s2
     }
