
posSsb =seq(length(res$names))[res$names=="ssb"]
posFbar=seq(length(res$names))[res$names=="fbar"]

cv=array(c(res$cov[rev(posFbar)[1],rev(posFbar)[1]],
           res$cov[rev(posSsb )[1],rev(posFbar)[1]],
           res$cov[rev(posSsb )[1],rev(posFbar)[1]],
           res$cov[rev(posSsb )[1],rev(posSsb)[1]]),dim=c(2,2), dimnames=list(c("fbar","ssb"),c("fbar","ssb")))

egn=eigen(cv)
inv=solve(egn$vectors)

circle=function(radius=1,interval=0.01) 
{ x=seq(0,1,interval)*radius
  x=c(-rev(x),x[-1])
  y=(radius*radius-x*x)^0.5
  
  #return(data.frame(x=x,y=y))
  x=c(x, rev(x[-1]), x[1])  
  y=c(y,-rev(y[-1]), y[1])
  
  return(as.matrix(data.frame(x=x,y=y)))}
cl=circle(qnorm(.75))
dat=as.data.frame(circle()%*%solve(cv))

ggplot(dat)+geom_path(aes(fbar,ssb))

cv=array(c(res$cor[rev(posFbar)[1],rev(posFbar)[1]],
           res$cor[rev(posSsb )[1],rev(posFbar)[1]],
           res$cor[rev(posSsb )[1],rev(posFbar)[1]],
           res$cor[rev(posSsb )[1],rev(posSsb)[1]]),dim=c(2,2), dimnames=list(c("fbar","ssb"),c("fbar","ssb")))


cv=        res$cov[c(rev(posFbar)[1],rev(posSsb)[1]),
                   c(rev(posFbar)[1],rev(posSsb)[1])]


STD = 2                     #%# 2 standard deviations
conf = pnorm(2)-pnorm(-2)   #%# covers around 95% of population
scale = chi2inv(conf,2)     #%# inverse chi-squared with dof=#dimensions
scale=dinvchisq(conf,2)
eigen(cv*scale)

Cov = cov(X0) * scale;
[V D] = eig(Cov);



corInv=solve(cv)


library(ellipse)
# Plot an approximate 95% confidence region for the Puromycin
# parameters Vm and K, and overlay the ellipsoidal region
data(Puromycin)
Purboth <- nls(formula = rate ~ ((Vm + delV * (state == "treated"))
                                 * conc)/(K + conc), data = Puromycin,
               start = list(Vm = 160, delV = 40, K = 0.05))
Pur.prof <- profile(Purboth)
plot(ellipse(Pur.prof, which = c('Vm', 'K')), type = 'l')
lines(ellipse(Purboth, which = c('Vm', 'K')), lty = 2)
params <- Purboth$m$getPars()
points(params['Vm'],params['K'])

fn=function(radius=.5,interval=0.0005,chiSq=2.3,corInv) 
{ 
  x=seq(0,1,interval)*.510
  y=((chiSq-x*x*(corInv[1,1]*corInv[2,2]+corInv[1,1]*corInv[1,1]))/(corInv[1,1]*corInv[2,2]+corInv[2,2]*corInv[2,2]))^0.5
  
  res=data.frame(x=x,y=y)
  
  res=rbind.fill(res,transform(res[rev(dimnames(res)[[1]]),],y=-y))
  
  res=rbind.fill(res,transform(res[rev(dimnames(res)[[1]]),],x=-x))
  
  ggplot(res)+geom_path(aes(x,y))
  
  
}
cl=circle(qnorm(.75))
dat=as.data.frame(circle()%*%solve(cv))
