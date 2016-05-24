"rrisk" <-
function(bound,pparms,sig.level=0.05,type="b"){
    d<-bound
    nd<-length(d$N)
    if (type=="b"){
        if (is.vector(pparms)){
            if (length(pparms)!=2){ stop("pparms is vector with length not equal to 2") }
            pparms<-matrix(pparms,1,2)
        }
        nout<-dim(pparms)[1]
        rr<-check<-EN<-rep(NA,nout)
        for (i in 1:nout){
            a<-pparms[i,1]
            b<-pparms[i,2]
            phi<-rep(0,nd)
            phi[d$p.value<=sig.level]<-1
            rr[i]<-sum( d$Kstar*beta.ratio(d$S+a,d$N-d$S+b,d$S+1,d$N-d$S+1)*(1/beta(a,b))*
                (phi+pbeta(sig.level,d$S+a,d$N-d$S+b)*(1-2*phi)) )
            check[i]<- sum( d$Kstar*beta.ratio(d$S+a,d$N-d$S+b,d$S+1,d$N-d$S+1)*(1/beta(a,b)) )
            EN[i]<- sum( d$N*d$Kstar*beta.ratio(d$S+a,d$N-d$S+b,d$S+1,d$N-d$S+1)*(1/beta(a,b)) )
         }
    }
    else if (type=="p"){
        nout<-length(pparms)
        rr<-check<-EN<-rep(NA,nout)
        for (i in 1:nout){
            if (pparms[i]<=sig.level) index<-(d$p.value>sig.level)
            else index<-(d$p.value<=sig.level)
            rr[i]<- sum(d$Kstar[index]*mult(d$N[index],d$S[index],pparms[i]) )
            check[i]<- sum( d$Kstar*mult(d$N,d$S,pparms[i]) )
            EN[i]<-sum(d$N*d$Kstar*mult(d$N,d$S,pparms[i]) )

        }
    }
    out<-list(check=check,rr=rr,EN=EN) 
  
    out
}

