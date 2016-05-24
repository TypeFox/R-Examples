`Boot.ts` <-
function(obj, R=1, ...){
p<-SelectModel(obj, lag.max=ceiling(length(obj)/4), Best=1)
ans<-FitAR(obj, 1:p)
Boot.FitAR(ans, R=R)
}

