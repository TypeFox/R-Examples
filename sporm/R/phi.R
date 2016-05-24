phi <-
function(N, theta, lambda) {
    alpha<-lambda-(1-lambda)*theta
    bta<-4*theta*lambda*(1-lambda)
    Ktheta<-sqrt(((theta-1)*((0:N)/N)+alpha)^2+bta)
    (1+((theta-1)*((2*(1:N)-1)/N)+2*alpha)/(Ktheta[-1]+Ktheta[-(N+1)]))/(2*lambda)
}

