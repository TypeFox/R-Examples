momeff.zeromean=function(tstat,n1,n2=n1,gamma2,lower.df=6.1,upper.df=100,approx=TRUE) {
########### mean zero effect size case:
    if(missing(n1) && !missing(gamma2) && gamma2>0){

    }else if (missing(n1)){ ## all, df and gamma2 need to be estimated, assuming n1=n2
        fobj.v=function(v){
            if(v<=6)return(NA)
            this.fit=parncpt.momeff(tstat=tstat,n1=(v+2)/2,zeromean=TRUE)
            var.ncp=attr(this.fit,'gamma2')*(v+2)/4
            (
            mean(tstat^6)-
            15*v^3/(v-2)/(v-4)/(v-6)*(this.fit+(1-this.fit)*(var.ncp^3+var.ncp^2+3*var.ncp+1))
            )^2
         }
         v=try(optimize(fobj.v,c(lower.df,upper.df)))
         if(inherits(v,'try-error'))return(NA)
         last.fit=parncpt.momeff(tstat,(v$minimum+2)/2,zeromean=TRUE)
         attr(last.fit,'df')=v$minimum
         attr(last.fit,'n1=n2')=(v$minimum+2)/2
         return(last.fit)
         
    }else if(missing(gamma2) || gamma2<=0){
        if(n1+n2<=6)return(NA)
        et2=mean(tstat^2);et4=mean(tstat^4);v=n1+n2-2
        pi0.num=(et4-3*et2^2)*v^2+(12*et2^2-6*et4)*v+8*et4-12*et2^2
        pi0.den=(et4-6*et2+3)*v^2+(12*et2-6*et4)*v+8*et4
        var.delta.den=(3*et2-3)*v*v-6*et2*v
        
        ans=pi0.num/pi0.den
        gamma2=pi0.den/var.delta.den*(n1+n2)/n1/n2
        attr(ans,'method')='zero mean, exact'
        if(ans<0 || ans>1){
            if(!approx)return(NA)
            boot.ans=boot.gam2=numeric(1000)
            for(r in 1:1000){## try bootstrap
                tstatb=sample(tstat,replace=TRUE)
                et2=mean(tstatb^2);et4=mean(tstatb^4);v=n1+n2-2
                pi0.num=(et4-3*et2^2)*v^2+(12*et2^2-6*et4)*v+8*et4-12*et2^2
                pi0.den=(et4-6*et2+3)*v^2+(12*et2-6*et4)*v+8*et4
                var.delta.den=(3*et2-3)*v*v-6*et2*v
                
                boot.ans[r]=pi0.num/pi0.den
                boot.gam2[r]=pi0.den/var.delta.den*(n1+n2)/n1/n2
            }
            
            if (all(boot.ans<0 | boot.ans>1)) { ## none of the bt sample works
                ans=if(ans<0) 0 else 1
            }else {                             ## some bt samples works
                pi0ok.idx=boot.ans>=0 & boot.ans<=1
                gam2ok.idx=boot.gam2>0
                ok.idx=pi0ok.idx & gam2ok.idx
                if(!any(ok.idx))ok.idx=pi0ok.idx
                ans=mean(boot.ans[ok.idx])
                gamma2=mean(boot.gam2[ok.idx])
            }
                
            attr(ans,'method')='zero mean, bootstrap approximate'
        }
        if(gamma2<0)warning('gamma2<0')
        
        attr(ans,'eff.mean')=0
        attr(ans,'gamma2')=gamma2

        return(ans)

    }else{## known gamma2
        if(n1+n2<=4)return(NA)
        
        ans=1-(mean(tstat^2)/(n1+n2-2)*(n1+n2-4)-1)/n1/n2/gamma2*(n1+n2)
        attr(ans,'method')='zero mean, exact, known gamma2'

        if(ans<0 || ans>1){
            if(!approx)return(NA)
            fobj=function(pi0)(
                mean(tstat^2)-(n1+n2-2)/(n1+n2-4)*(1+(1-pi0)*n1*n2*gamma2/(n1+n2))
                              )^2
            opt.fit=try(optimize(fobj,interval=c(0,1)))
            if(inherits(opt.fit,'try-error'))return(NA)
            ans=opt.fit$minimum
            attr(ans,'method')='zero mean, approximate, known gamma2'
        }
        attr(ans,'eff.mean')=0
        attr(ans,'gamma2')=gamma2
        return(ans)
    }
}
