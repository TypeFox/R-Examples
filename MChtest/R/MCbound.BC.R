"MCbound.BC" <-
function(Nmax,Smax,conf.level=.99){
    S<-c(rep(Smax,Nmax-Smax+1),(Smax-1):0)
    n<-c(Smax:Nmax,rep(Nmax,Smax))
    q<-Nmax-Smax+1
    Q<-length(S)
    p.value<- c(S[1:q]/n[1:q],(S[(q+1):Q]+1)/(n[(q+1):Q]+1) )
    alpha<-1-conf.level
    zero.one<-c(rep(1,q),rep(0,Q-q))
    #count<- choose(n-zero.one,S-zero.one)
    #cb<-count*beta(S+1,n-S+1)
    #check<- sum( count*(p^S)*((1-p)^(n-S) ) )
    cb<- c(S[1:q],n[(q+1):Q])/(n*(n+1))
    #check<-  sum( cb*mult(n,S,p=.03828483) )
    #print(check)
    ci.lower<-ci.upper<-rep(0,Q)
    alpha<- (1-conf.level)/2
    phat<-S/n
    eps<-.Machine$double.eps*2
    for (i in 1:Q){
        
        func1<-function(p,I=i){
            index<- phat >= phat[I]
            np<-length(p)
            out<-rep(NA,np)
            for (j in 1:np){
                 out[j]<-sum( cb[index]*mult(n[index],S[index],p[j]) )-alpha
            }
            out
        }
        if (i==Q){ ci.lower[i]<-0 }
        else { ci.lower[i]<-uniroot(func1,c(eps,1-eps))$root }
        func2<-function(p,I=i){
            index<- phat <= phat[I]
            np<-length(p)
            out<-rep(NA,np)
            for (j in 1:np){
                 out[j]<-sum( cb[index]*mult(n[index],S[index],p[j]) )-alpha
            }
            out
        }
        if (i==1) ci.upper[1]<-1
        else { ci.upper[i]<-uniroot(func2,c(eps,1-eps))$root }
    }
    parms<-c(Nmax,Smax)
    names(parms)<-c("Nmax","Smax")
    out<-list(S=S,N=n,p.value=p.value,ci.lower=ci.lower,
       ci.upper=ci.upper,Kstar=cb,conf.level=conf.level,type="BC",
       name="Sequential design (Besag and Clifford, 1991)",
       parms=parms)
    class(out)<-"MCbound"
    out
}

