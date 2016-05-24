`predict.FitAR` <-
function(object, n.ahead = 1, newdata = numeric(0),  ...){
    z<-object$z
    n<-length(z) #forecast origin
    z<-c(z,newdata)
    nr <- length(z)+n.ahead-1
    zm<-object$muHat
    phiHat<-object$phiHat
    sigsqHat<-object$sigsqHat
    r<-tacvfARMA(phi = phiHat, maxLag = nr, sigma2 = sigsqHat)
    TrenchForecast(z,r,zm,n,maxLead=n.ahead)
}

