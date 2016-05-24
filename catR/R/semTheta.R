semTheta<-function (thEst, it, x = NULL, model = NULL, D = 1, method = "BM", 
    priorDist = "norm", priorPar = c(0, 1), parInt = c(-4, 4, 
        33), constantPatt=NULL) 
{
constantPattern<-function(t) ifelse(sum(t)==0 | sum(t)==length(t),TRUE,FALSE)
METHOD<-NULL
if (!constantPattern(x) | !is.null(model) | is.null(constantPatt)) METHOD<-method
else{
if (sum(constantPatt==c("BM","EAP","WL"))==1) METHOD<-constantPatt
else RES<-Inf
}
if (!is.null(x)){
ind<-which(!is.na(x))
it<-it[ind,]
x<-x[ind]
}
if (!is.null(METHOD)){
    if (method == "EAP") {
        RES <- eapSem(thEst, it, x = x, model = model, D = D, 
            priorDist = priorDist, priorPar = priorPar, lower = parInt[1], 
            upper = parInt[2], nqp = parInt[3])
    }
    else {
if (is.null(model)) {
            dr0 <- function(th, it, D = 1, method = "BM", priorDist = "norm", 
                priorPar = c(0, 1)) {
                if (method == "BM") 
                  res <- switch(priorDist, norm = -1/priorPar[2]^2, 
                    unif = 0, Jeffreys = (sum(Ii(th, it, D = D)$d2Ii) * 
                      sum(Ii(th, it, D = D)$Ii) - sum(Ii(th, 
                      it, D = D)$dIi)^2)/(2 * sum(Ii(th, it, 
                      D = D)$Ii)^2))
                else res <- switch(method, ML = 0, WL = 0)
                return(res)
            }
            info <- sum(Ii(thEst, it, D = D)$Ii)
            res <- -dr0(thEst, it, D = D, method = method, priorDist = priorDist, 
                priorPar = priorPar) + info
            RES <- 1/sqrt(res)
        }
        else {
            met <- switch(method, ML = 1, BM = 2, WL = 3, EAP = 4)
            pd <- switch(priorDist, norm = 1, unif = 2, Jeffreys = 3)
            if (met == 1 | (met == 2 & pd == 2)) 
                optI <- sum(Ii(thEst, it, model = model, D = D)$Ii)
            if (met == 2 & pd == 1) 
                optI <- sum(Ii(thEst, it, model = model, D = D)$Ii) + 
                  1/priorPar[2]^2
            if ((met == 2 & pd == 3) | met == 3) {
                prI <- Ii(thEst, it, model = model, D = D)
                prJ <- Ji(thEst, it, model = model, D = D)
                if (met == 2) 
                  optI <- sum(prI$Ii) + (sum(prI$dIi)^2 - sum(prI$d2Ii) * 
                    sum(prI$Ii))/(2 * sum(prI$Ii)^2)
                else optI <- sum(prI$Ii)
            }
            RES <- 1/sqrt(optI)
        }
    }
}
    return(RES)
}
