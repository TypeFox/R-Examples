permcoefs.plsRnp <- function(dataRepYtt,ind,nt,modele,maxcoefvalues,wwetoile,ifbootfail){
dataRepYb=dataRepYtt[ind,1]
Tb=dataRepYtt[,-1]
tempCb=try(solve(t(Tb)%*%Tb)%*%t(Tb)%*%dataRepYb,silent=TRUE)
tempcoefs <- rbind(Intercept=0,wwetoile%*%tempCb)
    Cond <- FALSE
    try(Cond<-is.numeric(tempcoefs)&all(abs(tempcoefs)<maxcoefvalues),silent=TRUE)
    if (Cond) {
        return(tempcoefs)
    }
    else {
        return(ifbootfail)
    }
}
