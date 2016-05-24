tlmec <-
function(cens=NULL,y=NULL,x=NULL,z=NULL,nj=NULL,nu=4,family="t",criteria=TRUE,diagnostic=FALSE,initial,iter.max=200,error=1e-3){

require(mvtnorm)

if(length(y)==0) stop("the response is not supplied")
if(is.matrix(y)) y <- y[as.vector(!is.na(as.vector(t(y))))]
if(length(y) != sum(nj)) stop("not compatible sizes between the response y and the repetited measures nj")
if(length(cens)!= length(y)) stop("not compatible sizes between the response y and the censures cens")
if((family!="Normal") & (family!="t")) stop("family not recognized")

m<-length(nj)
initial<-inits(y,x,initial,ncol(as.matrix(z)),m)


if(family=="t") out <- EMT(cens,y,x,z,nj,nu,initial,criteria,diagnostic,error,iter.max)
if(family=="Normal") out <- EMN(cens,y,x,z,nj,initial,criteria,diagnostic,error,iter.max)

out
}

