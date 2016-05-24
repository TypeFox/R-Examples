pi.greco.grad <-
function(phi,i,ph.y,x.w,phi.tutte) 
{
phi.tutte[i,]<-phi
numobs<-nrow(x.w)
q.w<-ncol(x.w)
out<-colSums(ph.y[i,]*x.w)-(exp(phi%*%t(x.w))/colSums(exp(phi.tutte%*%t(x.w))))%*%x.w
return(out/numobs)
}
