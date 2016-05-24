coefs.plsRglmnp <- function(dataRepYtt, ind, nt, modele, family = NULL, maxcoefvalues,wwetoile, ifbootfail) 
{
dataRepYb=dataRepYtt[ind,1]
Tb=dataRepYtt[ind,-1]
tempCb=try(solve(t(Tb)%*%Tb)%*%t(Tb)%*%dataRepYb,silent=TRUE)
tempcoefs <- wwetoile%*%tempCb
    Cond <- FALSE
    try(Cond<-is.numeric(tempcoefs)&all(abs(tempcoefs)<maxcoefvalues),silent=TRUE)
    if (Cond) {
        return(tempcoefs)
    }
    else {
        return(ifbootfail)
    }
}
