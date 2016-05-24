momeff.nonzeromean=function(tstat,n1,n2=n1,gamma2,lower.df=6.1,upper.df=100,approx=TRUE) {
#################### non-zero mean effect size case

### some pre-computed constants
v=n1+n2-2
vv=v*v
et=mean(tstat)
et2=mean(tstat^2)
et3=mean(tstat^3)
etet=et*et
et2et2=et2^2

gv3=gamma(v/2-3/2)
sgv3=sqrt(gv3)
gv1=gamma(v/2-1/2)
gv=gamma(v/2)
gvgv=gv*gv
gv1gv1=gv1*gv1


#### after eliminating mu and sigma2, need to solve the following quadratic function of pi0
a=(
 -3*et*gv3*gv1gv1*v
 +2*et3*gv1gv1*gv1
)
b=(
 +3*et*et2*gv3*gv1gv1*v
 +3*et*gv3*gv1gv1*v
 -4*et3*gv1gv1*gv1
 -6*et*et2*gv3*gv1gv1
)
cc=(
 +2*et3*gv1gv1*gv1
 +6*et*et2*gv3*gv1gv1
 +4*et^3*gv3*gvgv
 -3*et*et2*gv3*gv1gv1*v
)
discriminant=b*b-4*a*cc

ans=ncp.mean=ncp.var=numeric(0)
if(discriminant>=0){ ## having real solutions
    pi01=(sqrt(discriminant)-b)/2/a ### solution 1
    mu=et/(1-pi01)/sqrt(v/2)/gv1*gv
    (sigma22=(et2-pi01*v/(v-2))/(1-pi01)/v*(v-2)-1-mu*mu)## this is based on 2nd moment eqn.
    (sigma23=(et3/(1-pi01)*2/v/sqrt(v/2)*gv/gv3-mu^3)/3/mu-1)## this is based on 3rd moment eqn. ## agrees:)

    if(pi01>=0 && pi01<=1 && sigma22>0 && abs(sigma22-sigma23)<1e-6){
        ans=c(ans,pi01)
        ncp.mean=c(ncp.mean,mu)
        ncp.var=c(ncp.var,(sigma22+sigma23)/2)
    }

    pi01=(-sqrt(discriminant)-b)/2/a ### solution 2
    mu=et/(1-pi01)/sqrt(v/2)/gv1*gv
    (sigma22=(et2-pi01*v/(v-2))/(1-pi01)/v*(v-2)-1-mu*mu)## this is based on 2nd moment eqn.
    (sigma23=(et3/(1-pi01)*2/v/sqrt(v/2)*gv/gv3-mu^3)/3/mu-1)## this is based on 3rd moment eqn. ## agrees:)

    if(pi01>=0 && pi01<=1 && sigma22>0 && abs(sigma22-sigma23)<1e-6){
        ans=c(ans,pi01)
        ncp.mean=c(ncp.mean,mu)
        ncp.var=c(ncp.var,(sigma22+sigma23)/2)
    }

    if(length(ans)==1){ ## having 1 feasible solution
        eff.mean=ncp.mean*sqrt(1/n1+1/n2)
        gamma2=ncp.var*(1/n1+1/n2)
        attr(ans,'eff.mean')=eff.mean
        attr(ans,'gamma2')=gamma2
        attr(ans,'method')='non-zero mean, exact solution'
        return(ans)
#        cat('unique pi0 ans=',ans,ncp.mean,sigma22,sigma23,fill=TRUE)
    }else if (length(ans)==2){ ## having 2 feasible solutions 
#        cat('double pi0 ans=',ans,ncp.mean,ncp.var,fill=TRUE)
        quants=quantile(tstat,seq(0,1,length=500))
        deltas=seq(quantile(tstat,.001),quantile(tstat,.999),length=201)
        dd=suppressWarnings(outer(quants,(deltas[-1]+deltas[-length(deltas)])/2,
                function(xx,yy)dt(xx,v,yy)))
        ll=c(-sum(log(
                ans[1]*dt(quants,v)+(1-ans[1])*dd%*%(pnorm(deltas[-1],ncp.mean[1],sqrt(ncp.var[1]))-
                                        pnorm(deltas[-length(deltas)],ncp.mean[1],sqrt(ncp.var[1])))
                )),-sum(log(
                ans[2]*dt(quants,v)+(1-ans[2])*dd%*%(pnorm(deltas[-1],ncp.mean[2],sqrt(ncp.var[2]))-
                                        pnorm(deltas[-length(deltas)],ncp.mean[2],sqrt(ncp.var[2])))
                )))
        ans=ans[which.min(ll)]
        eff.mean=ncp.mean[which.min(ll)]*sqrt(1/n1+1/n2)
        gamma2=ncp.var[which.min(ll)]*(1/n1+1/n2)
        attr(ans,'eff.mean')=eff.mean
        attr(ans,'gamma2')=gamma2
        attr(ans,'method')='non-zero mean, exact solution'
        return(ans)
#        cat('mle pi0 ans=',ans,ncp.mean[which.min(ll)],ncp.var[which.min(ll)],fill=TRUE)
    }
}

## if program runs thru here, neither solutions are feasible
if(!approx)return(NA)

    opt.fit=try(optimize(function(pp)abs(a*pp*pp+b*pp+cc),
                c(0,lastbin(2*pt(abs(tstat),v,lower.tail=FALSE))), ### this range could be changed to c(0,1)
                tol=1e-12))
    if(inherits(opt.fit,'try-error')) return(NA)
    
#    print(opt.fit)
    pi01=opt.fit$minimum
    mu=et/(1-pi01)/sqrt(v/2)/gv1*gv
    (sigma22=(et2-pi01*v/(v-2))/(1-pi01)/v*(v-2)-1-mu*mu)## this is based on 2nd moment eqn.
    (sigma23=(et3/(1-pi01)*2/v/sqrt(v/2)*gv/gv3-mu^3)/3/mu-1)## this is based on 3rd moment eqn. 
    ## BECAUSE THIS IS APPROX SOLUTION, THE ABOVE TWO WILL NOT AGREE. 
#    cat(pi01,mu,sigma22,sigma23,fill=TRUE)
    if(sigma22>0 && sigma23>0 ){## two solutions, check mle
        ans=rep(pi01,2)
        ncp.mean=rep(mu,2)
        ncp.var=c(sigma22,sigma23)
#        cat('double approx pi0 ans=',ans,ncp.mean,ncp.var,fill=TRUE)
        quants=quantile(tstat,seq(0,1,length=500))
        deltas=seq(quantile(tstat,.001),quantile(tstat,.999),length=201)
        dd=suppressWarnings(outer(quants,(deltas[-1]+deltas[-length(deltas)])/2,
                function(xx,yy)dt(xx,v,yy)))
        ll=c(-sum(log(
                ans[1]*dt(quants,v)+(1-ans[1])*dd%*%(pnorm(deltas[-1],ncp.mean[1],sqrt(ncp.var[1]))-
                                        pnorm(deltas[-length(deltas)],ncp.mean[1],sqrt(ncp.var[1])))
                )),-sum(log(
                ans[2]*dt(quants,v)+(1-ans[2])*dd%*%(pnorm(deltas[-1],ncp.mean[2],sqrt(ncp.var[2]))-
                                        pnorm(deltas[-length(deltas)],ncp.mean[2],sqrt(ncp.var[2])))
                )))
        ans=ans[which.min(ll)]
        eff.mean=ncp.mean[which.min(ll)]*sqrt(1/n1+1/n2)
        gamma2=ncp.var[which.min(ll)]*(1/n1+1/n2)
        attr(ans,'eff.mean')=eff.mean
        attr(ans,'gamma2')=gamma2
        attr(ans,'method')='non-zero mean, approximate solution'
        return(ans)
#            cat('mle approx pi0 ans=',ans[1],ncp.mean[1],ncp.var[which.min(ll)],fill=TRUE)
    }else if (sigma22*sigma23<0){ ## one solution
        ans=pi01
        ncp.mean=mu
        ncp.var=max(c(sigma22,sigma23))
        eff.mean=ncp.mean*sqrt(1/n1+1/n2)
        gamma2=ncp.var*(1/n1+1/n2)
        attr(ans,'eff.mean')=eff.mean
        attr(ans,'gamma2')=gamma2
        attr(ans,'method')='non-zero mean, approximate solution'
        return(ans)
#            cat('unique approx pi0 ans=',ans,ncp.mean,sigma22,sigma23,fill=TRUE)
    }else { ## no solution, report result, but giving a warning
        ans=pi01
        ncp.mean=mu
        ncp.var=mean(c(sigma22,sigma23))
        warning('gamma2<0')
        eff.mean=ncp.mean*sqrt(1/n1+1/n2)
        gamma2=ncp.var*(1/n1+1/n2)
        attr(ans,'eff.mean')=eff.mean
        attr(ans,'gamma2')=gamma2
        attr(ans,'method')='non-zero mean, approximate solution'
        return(ans)
#            cat('bad approx pi0 ans=',ans,ncp.mean,sigma22,sigma23,fill=TRUE)
   }

}## of non-zero mean model

