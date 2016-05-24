`chooseMatrix` <-
function(n,m){
    n.ch.m<-choose(n,m)
    if (n.ch.m==1){
        if (n==m) cm<-matrix(rep(1,n),1,n)
        else if (m==0) cm<-matrix(rep(0,n),1,n)
    }
    else{
        cm<-matrix(0,n.ch.m,n)
        cm[1,]<-c(rep(1,m),rep(0,n-m))
        next.row<-function(x,N=n){ 
            max.zero<-max((1:N)[x==0])
            max.one<-max((1:max.zero)[x[1:max.zero]==1])
            msum<-sum(x[(max.one):N])
            if (msum>1){
                x<-c(x[0:(max.one-1)],0,rep(1,msum),rep(0,n-msum-max.one))
            }
            else{
                x[max.one]<-0
                x[max.one+1]<-1
            }
            x
        }
        for (i in 2:n.ch.m){
            cm[i,]<-next.row(cm[i-1,])
        }
    }
    cm
}

