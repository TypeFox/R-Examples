cumsum.test<-function(x, alternative = 
        c("one.change", "epidemic"))
{
    Sn<-0
    if(is.vector(x)){
        n<-length(x)
        x<-matrix(x,n,1)
    }
    if(!is.matrix(x)) x<-as.matrix(x)
    if(is.matrix(x)) n<-nrow(x)
    Sigma<-var(x)
    alternative <- match.arg(alternative)
    ifelse(alternative == "epidemic", nc<-2, nc<-1)
    for(k in 1:(n-1))
        for(m in ((nc == 1)*n + (nc == 2)*(k+1)):n){
#            tmp<-abs(sum(x[(k+1):m]-(m-k)*sum(x)/n))
            tmp<-apply(as.matrix(x[(k+1):m,]),2,sum)-(m-k)*apply(x,2,sum)/n
            Tmp<-t(tmp)%*%solve(Sigma)%*%tmp
            if(nc==1){
                Tmp<-n*Tmp/((m-k)*(n-m+k))
                if(Sn<Tmp){
                    Sn<-Tmp
                    k.hat<-k
                }
            }
            else{
                Tmp<-Tmp/n
                if(Sn<Tmp){
                    Sn<-Tmp
                    k.hat<-k
                    m.hat<-m
                }
            }
        }
    Sn<-Sn^2/(n*var(x))
    if(nc==1) list(Sn = Sn, k.hat = k.hat)
    else list(Sn = Sn, k.hat = k.hat, m.hat = m.hat)
}
