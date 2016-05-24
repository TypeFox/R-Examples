SemiparChangePoint<-function(x, alternative = c("one.change", "epidemic"),
        adj.Wn = FALSE, tol=1.0e-7, maxit=50,trace=FALSE,... )
{
    nc<-0; n<-0
    alternative <- match.arg(alternative)
    ifelse(alternative == "epidemic", nc<-2, nc<-1)
    alpha.hat<-0;
    if(!is.matrix(x)) x<-as.matrix(x)
    if(is.matrix(x)){
        n<-nrow(x)
        r<-ncol(x)
    }
    if(is.vector(x)){
        n<-length(x); r<-1
        x<-matrix(x,n,1)
    }
    if(nc ==1 ){
        ll<-NULL
        ll[1]<--n*log(n); ll[n+1]<-ll[1]
    }
    else{
        ll<-matrix(0,n,n)
        ll[1,n]<--n*log(n)
        for(i in 1:n){
            ll[i,i]<--n*log(n);
        }
    }
    Beta<-NULL; z0<-NULL
    beta.hat<-NULL
    BETA<-NULL
    lik.temp<--n*log(n); k.hat=n; m.hat=n
    Vn<-0; Wn<-0;it<-0
    Sn<-0; ind<-NULL; cnt<-0
    for(k in 1:(n-1)){
        for(m in ((nc == 1)*n + (nc == 2)*(k+1)):n){
            it<-it+1
            k1<-k+1; m1<-m+1; nm<-n-m
            mk<-m-k;
            z0[1:k]<-rep(0,k)
            if(m<n) z0[m1:n]<-rep(0,nm)
            z0[k1:m]<-rep(1,mk)
            z.glm<-glm(formula=z0~x,family=binomial,control=glm.control(tol,maxit,trace))
            Beta<-coef(z.glm)
            rho<-(k+nm)/mk
            ALPHA<-Beta[1]+log(rho); BETA<-as.vector(Beta[2:(r+1)])
            temp1<-sum(log(mk*exp(ALPHA+x%*%BETA)+(k+nm)))
            temp2<-sum(ALPHA+as.matrix(x[k1:m,])%*%BETA)
            if(nc ==1){ 
                ll[k+1]<--temp1+temp2;
                if(ll[k+1]==-Inf){
                    cnt<-cnt+1
                    ind[cnt]=k+1
                }
                if((ll[k+1]>lik.temp)||(ll[k+1]>=lik.temp && k<k.hat)){
                    lik.temp<-ll[k+1]
                    k.hat<-k
                    alpha.hat<-ALPHA
                    beta.hat<-Beta[2:(r+1)]
                }
            }
            else{
                ll[k,m]<--temp1+temp2; ll[m,k]<-ll[k,m]
                if((ll[k,m]>lik.temp)||(ll[k,m]>=lik.temp && m-k>m.hat-k.hat)){
                    lik.temp<-ll[k,m]
                    k.hat<-k; m.hat<-m;
                    alpha.hat<-ALPHA
                    beta.hat<-Beta[2:(r+1)]
                }
                Lambda<-(ll[k,m]-ll[1,1])
                Vn<-Vn+2*(6*(n-m+k)*(m-k)-1)*Lambda/(3*n^4)
                if(Lambda==-Inf) Lambda<-0
                ctemp<-max((m-k)*(n-m+k)*Lambda,(m+1-k)*(n-m+k-1)*Lambda,(m-k-1)*(n-m+k+1)*Lambda)
                if(Wn<2*ctemp/(n*n)){
                    Wn<-2*ctemp/(n*n)
                    if(adj.Wn && r == 1) Wn <- Wn*(1+0.155/sqrt(n)+0.24/n)
    #               cat("Iteration ",it,"(k,m)=(",k,",",m,")\n")
                }
            }
        }
    }
    if(nc==1){
        Sn<-2*(lik.temp-ll[1])
        if(cnt>0) ll<-replace(ll,ind,rep(min(ll[ll!=-Inf]),cnt))
        return (list(k.hat = k.hat, ll = ll, Sn = Sn,  alpha.hat = alpha.hat, beta.hat = beta.hat))
    }
    if(nc==2) list(k.hat = k.hat ,m.hat = m.hat, ll = ll, Vn = Vn, Wn = Wn, alpha.hat = alpha.hat, beta.hat = beta.hat)    
}
