getbetap.A <-
function(getAmoment.list, A=NULL, fix.obs=F){

# this function is derived from the ideas in getbetap.new
 # It assumes each A has been forced to have mean zero
 x=getAmoment.list$x
 y=getAmoment.list$y
 n=getAmoment.list$n
 z=getAmoment.list$z
 if (length(A)==0){A=getAmoment.list$A}
 else {A=A-getAmoment.list$mu}
 EA=getAmoment.list$EA # this should be a vector of zeros
 EA2=getAmoment.list$EA2
 EA3=getAmoment.list$EA3
 EA4=getAmoment.list$EA4
 V=EA2-EA^2
 s=(EA3-3*EA*V-3*EA^3)/V^(3/2)
 k=EA4/V^2-3
 effective.n= -1-6/(-abs(k))
 lowest.alpha=2/factorial(effective.n) # this is an estimate of the smallest achievable p-value
 lowest.alpha[EA2==0]=1 # if there is no variance in A, then clearly the smallest p-value is 1

 r=A/sqrt(V*(n-1))
 alpha1= (3*k + 36*s*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) - 18*s^3*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) - 3*s^2 + 3*k^2*s*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) - 3*k*s^3*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) + 24*k*s*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) + 6)/(2*k - 3*s^2) - (- 6*s^2 + 6*k + 12)/(2*k - 3*s^2)
 alpha2=(3*k - 36*s*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) + 18*s^3*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) - 3*s^2 - 3*k^2*s*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) + 3*k*s^3*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) - 24*k*s*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) + 6)/(2*k - 3*s^2) - (- 6*s^2 + 6*k + 12)/(2*k - 3*s^2)
 beta1= -(3*k + 36*s*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) - 18*s^3*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) - 3*s^2 + 3*k^2*s*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) - 3*k*s^3*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) + 24*k*s*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) + 6)/(2*k - 3*s^2)
 beta2= -(3*k - 36*s*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) + 18*s^3*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) - 3*s^2 - 3*k^2*s*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) + 3*k*s^3*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) - 24*k*s*(-1/(- k^2*s^2 + 32*k^2 - 84*k*s^2 + 96*k + 36*s^4 - 180*s^2))^(1/2) + 6)/(2*k - 3*s^2)   
 alpha=alpha1
 beta=beta1
 which.negative=grep(T,(alpha<=0)|(beta<=0))
 alpha[which.negative]=alpha2[which.negative]
 beta[which.negative]=beta2[which.negative]
 beta.mean=alpha/(alpha+beta)
 beta.var=(alpha*beta)/((alpha+beta)^2*(alpha+beta+1))
 c0=beta.mean
 c1=sqrt(beta.var*(n-1))

 rprime=c0+c1*r

 rprime.high=c0+c1*abs(r)
 rprime.low=c0-c1*abs(r)

 twosidedp=pbeta(rprime.high,alpha,beta,lower.tail=F)+pbeta(rprime.low,alpha,beta)
 rightp=pbeta(rprime,alpha,beta,lower.tail=F)
 leftp=pbeta(rprime,alpha,beta)
 doublep=2*apply(cbind(rightp,leftp),1,min)

# compute Z
 t=A/sqrt(EA2)
 #pt=2*pt(-abs(t),df=effective.n-2)
 r2=t^2/(n-2+t^2)
 pt=pbeta(r2,1/2,.5*(n-2),lower.tail=F)

 if (fix.obs==T){
   which=unique(c(grep(T,rightp==0),grep(T,leftp==0),grep(T,twosidedp==0)))
   if (length(which)>1){
     x0=x[which,]
     result0=getbetap.A.2(x0,y,z=z)
     leftp[which]=result0$leftp
     rightp[which]=result0$rightp
     twosidedp[which]=result0$twosidedp
     doublep[which]=2*apply(cbind(rightp[which],leftp[which]),1,min)
     }
     }
 
# find p-values failing Chebyshev, causing them to be too large. Set them equal to Chebyshev bound
 chebyshev.p=EA2/A^2
 fail.chebyshev=grep(T,(chebyshev.p<twosidedp))
 twosidedp[fail.chebyshev]=chebyshev.p[fail.chebyshev]
 doublep[fail.chebyshev]=chebyshev.p[fail.chebyshev]
 leftp[fail.chebyshev]=chebyshev.p[fail.chebyshev]/2
 rightp[fail.chebyshev]=chebyshev.p[fail.chebyshev]/2

# clean up any remaining zero p-values using Spearman
 if (fix.obs==T){
   which=unique(c(grep(T,rightp==0),grep(T,leftp==0),grep(T,twosidedp==0)))
   if (length(which)==1){
     current=cor.test(x[which,],y,method="spearman")
     doublep[which]=twosidedp[which]=current$p.value
    }  
   if (length(which)>1){
     for (j in (1:length(which))){
     current=cor.test(x[which[j],],y,method="spearman")
     doublep[which[j]]=twosidedp[which[j]]=current$p.value
     }
     }
      }

 return(list(twosidedp=twosidedp,rightp=rightp,leftp=leftp,pdouble=doublep,chebyshev.p=chebyshev.p,pt=pt,lowest.alpha=lowest.alpha))
 }
                                       
