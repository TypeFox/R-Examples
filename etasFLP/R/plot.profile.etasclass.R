plot.profile.etasclass=function(x,prob=c(0.90,0.95,0.99),...){
## plot profile likelihood of an etas model
# 
if(sum(is.element(class(x),"profile.etasclass"))==0) stop("object is not of the required class profile.etasclass")

namespar=c("mu","k0","c","p","a","gamma","d","q")

xx=x$param.vec
yy=x$logl.vec
zz=2*(yy-min(yy))
namepar=namespar[x$iprofile]

perc=qchisq(prob,df=1)

sn=spline(xx,zz,n=200,method="natural")
plot(sn,type="l",xlab=namepar,main=paste("Profile -2[log(LR)] for parameter ",namepar),
	  ylab="-2[log(LR)]")
abline(h=perc,col="red")
grid(col="darkgrey")
text(x=xx[(1+length(xx))/2],y=perc,label=prob)

### Inverse spline interpolation for asymptotic confidence intervals
n	=length(xx)
nc	=trunc((n+1)/2)
x1	=xx[1:nc]
y1	=zz[1:nc]

sf=splinefun(y1,x1,method="natural")
l1=round(sf(perc),3)

x2=xx[nc:n]
y2=zz[nc:n]

sf=splinefun(y2,x2,method="natural")
l2=round(sf(perc),3)

conf=cbind(prob,l1,l2)
colnames(conf)=c("Coverage","Lower","Upper")
rownames(conf)=1:length(prob)

cat("Asymptotic  confidence  intervals:","\n")
print(conf)

ris=list(spline.profile=sn,conf=conf,prob=prob)
}


