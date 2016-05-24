V.theta <-
function(t, N, theta, lambda){
    alpha<-lambda-(1-lambda)*theta
    bta<-4*theta*lambda*(1-lambda)
    sqrt(((theta-1)*t+alpha)^2+bta)
}

