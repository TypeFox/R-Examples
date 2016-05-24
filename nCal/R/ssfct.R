# self start functions
ssfct.drc.1.5.2 <- function (dataFra) {
    x <- dataFra[, 1]
    y <- dataFra[, 2]
    startVal <- rep(0, 5)
    startVal[3] <- max(y) + 0.001
    startVal[2] <- min(y) - 0.001
    startVal[5] <- 1
    indexT2 <- (x > 0)
    x2 <- x[indexT2]
    y2 <- y[indexT2]
    startVal[c(1, 4)] <- find.be2(x2, y2, startVal[2] - 0.001, 
        startVal[3])
#    print (startVal)
    return(startVal)
}
find.be2 <- function(x, y, c, d)
{
#    myprint (x," ")
#    myprint (y)
#    myprint (c)
#    myprint (d)
    logitTrans <- log((d - y)/(y - c))  
    
    lmFit <- lm(logitTrans ~ log(x))
#        eVal <- exp((-coef(logitFit)[1]/coef(logitFit)[2]))
#        bVal <- coef(logitFit)[2]
    
    coefVec <- coef(lmFit)
    bVal <- coefVec[2]        
    eVal <- exp(-coefVec[1]/bVal)    
    
    return(as.vector(c(bVal, eVal)))
}
ss.fct.via.LL4= function (dataFra) {
    x <- dataFra[, 1]
    y <- dataFra[, 2]
    fit0= drm(y ~ x, fct = LL.4(), weights= y^-.5)
    start=c(coef(fit0), 1)
    return (start)
}
