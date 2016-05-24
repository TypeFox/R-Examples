"ahaz.adjust"<-function(surv,X,weights,idx,method=c("coef","z","crit"))
  {
    ## Purpose: Fast calculation of univariate measures of association for columns in X,
    ##          adjusted for covariates specified in 'idx'
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   surv       : Surv object (right censored/counting process)
    ##   X          : Numeric matrix of covariates
    ##   weights    : Weight vector (nonnegative)
    ##   idx        : Which columns to adjust for?
    ##   method     : Adjustment method
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen

     this.call <- match.call()
     method<-match.arg(method)
     
     if(missing(idx) || !is.numeric(idx) || any(idx<1) || any(idx>ncol(X)))
       stop("'idx' incorrectly specified")
     idx<-as.integer(idx)
     
     idx <- sort(idx)

     if(method == "z")
       {
         # Partial solutions for covariates in 'idx
         m<-ahaz.partial(surv, X, weights, idx = idx)
         Dbig<-m$D
         Bbig<-m$B

         # Covariates in 'idx' only
         Borig<-Bbig[,idx]
         Dorig<-Dbig[,idx]

         # Univariate results
         muni<-ahaz(surv,X, weights, univariate=TRUE)
         d<-muni$d
         D<-muni$D
       
         
         iA<-solve(Dorig)
         n<-ncol(iA)
         vf<-iA%*%Dbig
         ik<-1/(D-apply(Dbig*vf,2,sum))
         x1<--ik*t(vf)
         nomin<-x1%*%d[idx]+ik*d

         denomin<-ik^2*muni$B+2*ik*apply(Bbig*t(x1),2,sum)+apply(t(x1) * (Borig%*%t(x1)), 2, sum)
         denomin[idx]<-NA
         nomin[idx]<-NA
         adj<-nomin/sqrt(denomin)
          # return(adj)
       } else  {
         m1<-ahaz.partial(surv=surv,X=X,weights=weights, idx=idx)
         m3<-ahaz(surv=surv,X=X,weights=weights, univariate=TRUE)
         
         DD<-m1$D
         iA<-solve(DD[,idx])
         D<-m3$D
         d<-m3$d
         
         n<-ncol(iA)
         vf<-iA%*%DD
         ik<-1/(D-apply(DD*vf,2,sum))
         x1<--ik*t(vf)%*%d[idx]
         x2<-ik*d
         adj<-x1+x2
         adj[idx]<-NA
         if(method!="coef"){
           adj<--(ik*d^2+2*d*x1+x1^2/ik)
           adj[idx]<-NA
         } 
       }
     return(list("call"=this.call,"idx"=idx,"adj"=drop(adj)))
   }
