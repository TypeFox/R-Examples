pi.greco.hess <-
function(phi,i,ph.y,x.w,phi.tutte) 
{
phi.tutte[i,]<-phi
numobs<-nrow(x.w)
a<-(exp(phi%*%t(x.w))/colSums(exp(phi.tutte%*%t(x.w))))
out<-(t(x.w)%*%x.w)%o%(a^2-a)
out<-apply(out[,,1,],1,rowMeans)
return(out/numobs)
}
