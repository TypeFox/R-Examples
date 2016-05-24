ks.stat <-
function(x, y){
    m<-length(x)
    n<-length(y)
    N<-m+n
    lambda<-m/N
    z<-c(y,x)
    z<-sort(z)
    Fx<-ecdf(x)
    Gy<-ecdf(y)
    theta0<-dd.est(x,y)
    p0<-phi(N, theta0, lambda)/N
    res<-try(mrle.sporm(x, y, theta0, p0), TRUE)
    theta<-res$theta
    p<-res$p
    F.tilde<-cumsum(p)
    G.tilde<-theta*F.tilde/(1+(theta-1)*F.tilde)
    ks<-sqrt(m)*max(abs(F.tilde-Fx(z)))+sqrt(n)*max(abs(G.tilde-Gy(z)))
    list(ks=ks, theta=theta)
}
