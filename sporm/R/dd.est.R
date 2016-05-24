dd.est <-
function(x, y){
    r<-3
    m<-length(x)
    n<-length(y)
    Ni<-m*(1-ecdf(x)(sort(y)))
    Mi<-n+1-(1:n)
    theta0.hat<-(r-1)*n*sum(Ni^r/Mi^2)/m^r
    theta0.hat
}
