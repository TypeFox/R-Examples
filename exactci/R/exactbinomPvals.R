`exactbinomPvals` <-
function(x,n,p,relErr=1+10^(-7),
    tsmethod="minlike",midp=FALSE){

    np<-length(p)
    if (any(p>1) | any(p<0)) stop("p must be in [0,1]")
    #fisher<-blaker<-central<-rep(NA,nor)
    pvals<-rep(NA,np)
    support<- 0:n
    X<- support==x
    ns<-length(support)
    if (tsmethod=="blaker"){
        for (i in 1:np){
            f<-dbinom(support,n,p[i])
            d<-f[X]
            F<-cumsum(f)
            Fbar<-cumsum(f[ns:1])[ns:1]
            if (midp){
                F<- F - 0.5*f
                Fbar<- Fbar - 0.5*f
            }
            lower<-F[X]
            upper<-Fbar[X]
            if (lower<=upper){
                if (any(Fbar<=lower*relErr)){
                    pvals[i]<-lower + 
                      max(Fbar[Fbar<=lower*relErr])
                } else pvals[i]<-lower
            } else if (upper<lower){
                if (any(F<=upper*relErr)){
                    pvals[i]<-upper + 
                      max(F[F<=upper*relErr])
                } else pvals[i]<-upper
            } 
        }
    } else if (tsmethod=="minlike"){
        for (i in 1:np){
            f<-dbinom(support,n,p[i])
            d<-f[X]
            if (midp){
                pvals[i]<-sum(f[f<=d*relErr])-0.5*d
            } else {
                pvals[i]<-sum(f[f<=d*relErr])
            }
        }   
    } else if (tsmethod=="central"){
        Xlo<-support<=x
        Xhi<-support>=x
        for (i in 1:np){
            f<-dbinom(support,n,p[i])
            if (midp){
                pvals[i]<-2*min(sum(f[Xlo])-0.5*f[X],
                                sum(f[Xhi])-0.5*f[X])
            } else {
                pvals[i]<-2*min(sum(f[Xlo]),sum(f[Xhi]))
            }
        }
    }
    pvals<-pmin(1,pvals)
    out<-list(p=p,pvals=pvals)
    out
}
