ELF <-
function(xdat, grp, nperm=50, normalize=FALSE)
{
 #require(rsmooth)  # download rsmooth from www.meb.ki.se/~yudpaw

 if (nperm==0) {normalize=TRUE}
 ng = nrow(xdat)
 stat = tstatistics2(xdat, grp)$tstat
 pval = 2*(1-pt(abs(stat), df = ncol(xdat)-2))
 if (normalize) stat = qnorm(pt(stat, df = ncol(xdat)-2))

# patterns from SVD
 fun = ufun(xdat, grp, nperm=nperm, normalize=normalize)
 delta = fun$delta
 zk = fun$zk
 zz = c(zk[1]-delta/2, zk+ delta/2)
 phi0 = ng*delta*fun$phi0
 phi1 = ng*delta*fun$u1

# first step
 y = table(cut(stat,zz))
 yreg = lm(y~ -1+phi0 + phi1)
# second step
 res = y - yreg$fit
 #f1 = rsmooth(zk,res, lam=1)$y
f1 = lowess(zk,res)$y
   f1 = ifelse(f1<0,0,f1)
   yreg = lm(y~ -1+phi0 + phi1 + f1)

# parameters
 p0.elf = yreg$coef[1]
 b = yreg$coef[2]/yreg$coef[1]

# fitted null 
 sig = 1 + b/sqrt(2)
 y0= ng*delta* dnorm(zk, sd = sig)
 #if (abs(b)< 0.2) y0= ng*delta* dnorm(zk, sd = sig)
 #if (abs(b)> 0.2) y0 =  phi0 + b*phi1
 
 #if (min(y0)<0) y0= ng*delta* dnorm(zk, sd = sig)
 #if (min(y0)<0) y0= pmax(y0,0)
 #y0= ng*delta*dt(zk/sig, df=10)/sig  
 
 if (nperm>0 && abs(b)>0.2) {fun = ufun(xdat, grp, nperm=nperm, normalize=FALSE);lim=0.5;p0.elf<- min(sum(pval>lim)/((1-lim)*ng),1);
 err=fun$Y-as.vector(y);sumerr= apply(err^2*as.vector(y),2,sum);closest = order(sumerr)[1:5];y0 = apply(fun$Y[,closest],1,mean)}
  

# ................................. FDR computation: right-hand side
rej = rev(cumsum(rev(y[zk>0])))
  #rej = ifelse(rej<1, 1, rej)  # avoid zero
  rej0 = p0.elf*rev(cumsum(rev(y0[zk>0])))  # rejecti
elf = pmin(rej0/ rej,1)
  elfp = cummin(elf)   # monotone constrain

# ................................. left-hand side
rej = cumsum(y[zk<= 0])
  #rej = ifelse(rej<1, 1, rej)  # avoid zero
  rej0 = p0.elf*cumsum(y0[zk<= 0])  # null rejection
elf = pmin(rej0/ rej,1)
  elfm = rev(cummin(rev(elf)))   # monotone constrain

# combined
elf = c(elfm,elfp)

# ELF for individual genes
indiv.elf = approx(zk,elf, xout=stat, rule=2)$y

return(list(zk=zk,y=y, funELF=elf, 
       stat=stat, ELF= indiv.elf, p0= p0.elf, 
       b=b, phi0=phi0, phi1=phi1, delta=delta, y0=y0))
}
