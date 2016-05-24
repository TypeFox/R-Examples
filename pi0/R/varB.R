varB=function(m,deltamax,K){
    stopifnot(m>=1 && m<=K-1)
    (6*m*m-12*m+7)/6*(deltamax/(K-1))^2
}
