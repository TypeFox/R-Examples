`effectsize` <-
function(tstat,ntilde,m)
{
# cm=gamma(m/2)/(sqrt(m/2)*gamma((m-1)/2))
ln.cm<-lgamma(m/2)-log(sqrt(m/2))-lgamma((m-1)/2)  
cm<-exp(ln.cm)
d=tstat/sqrt(ntilde)
dprime=cm*d
terme1=m/((m-2)*ntilde)
vard=terme1+d^2*(terme1*ntilde-1/cm^2)
vardprime=cm^2*(terme1+dprime^2*(terme1*ntilde-1/cm^2))
result=cbind(d,vard,dprime,vardprime)
colnames(result)=c("d","vard","dprime","vardprime")
result
}

