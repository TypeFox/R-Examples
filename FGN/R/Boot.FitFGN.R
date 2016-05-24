`Boot.FitFGN` <-
function(obj, R=1, ...){
n<-length(obj$res)
H<-obj$H
sig<-sqrt(obj$sigsq)
zm<-obj$muHat
if (R==1){
    z<-zm+SimulateFGN(n, H)*sig
    if (is.null(obj$tsp))
        z
    else
        ts(z, start=obj$tsp[1], frequency=obj$tsp[3])
    }
else {
     x<-matrix(0, nrow=n, ncol=R)
     for (i in 1:R)
        x[,i]<-zm+SimulateFGN(n, H)*sig
    x
    }
}

