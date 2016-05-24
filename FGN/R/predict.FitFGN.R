`predict.FitFGN` <-
function(object, n.ahead = 1, ...){
    z<-object$z
    n<-length(z)
    H<-object$H
    zm<-object$muHat
    r<-var(z)*acvfFGN(H, n+n.ahead-1)
    TrenchForecast(z,r,zm,n,maxLead=n.ahead)
}

