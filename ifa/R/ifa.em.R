"ifa.em" <-function(y,ni,it=15,eps=0.001,init=NULL,scaling=TRUE)
{
if (scaling) y<-scale(y)
ni<-matrix(ni)
L<-nrow(ni)
numobs<-nrow(y)
numvar<-ncol(y)
ybar<- apply(y, 2., mean)
y<-scale(y, ybar, scale=FALSE) 
totni<-prod(ni)
maxni<-max(ni)

if (is.null(init)) output.init<-ifa.init.pca(y,L) else output.init<-init

psi<-output.init$psi
psi<-diag(diag(psi))
H<-output.init$H


w<-matrix(0,L,maxni)
mu<-matrix(0,L,maxni)
vu<-matrix(0,L,maxni)

for (i in 1:L) for (j in 1:ni[i]) {w[i,j]<-runif(1,0,1)
                                        mu[i,j]<-runif(1,-1,1)
                                        vu[i,j]<-runif(1,0,1)} 
w<-w/rowSums(w)

likelihood<-matrix(0,it)

storage.mode(L)<-"integer"
storage.mode(numobs)<- "integer"
storage.mode(numvar)<- "integer"
storage.mode(totni)<- "integer"
storage.mode(maxni)<- "integer"
storage.mode(ni) <- "integer"
storage.mode(it) <- "integer"

storage.mode(H) <- "double"
storage.mode(psi) <- "double"
storage.mode(mu) <- "double"
storage.mode(vu) <- "double"
storage.mode(y) <- "double"
storage.mode(w) <- "double"
storage.mode(eps) <-"double"
storage.mode(likelihood)<-"double"
pqy<-matrix(0,numobs,totni)
storage.mode(pqy)<-"double"
sigma<-array(0,c(L,L,totni))
storage.mode(sigma) <- "double"
EExxy<-matrix(0,L,L)
storage.mode(EExxy) <- "double"


out<-.Fortran("ifaem",y,numobs,L,numvar,ni,totni,maxni,it,H,w,mu,vu,eps,psi,likelihood,sigma,pqy,EExxy,PACKAGE="ifa")

likelihood<-out[[15]]
sigma<-out[[16]]
pqy<-out[[17]]
H<-out[[9]]
w<-out[[10]]
mu<-out[[11]]
vu<-out[[12]]
psi<-out[[14]]
psi<-diag(diag(psi))
likelihood<-matrix(likelihood[!likelihood==0])
EExxy<-out[[18]]


niter<-length(likelihood)
temp<-(likelihood-likelihood[niter])
r<-log(temp[-1]/temp[-niter])
r<-exp(mean(r[1:(niter-2)],na.rm=T))

se.H<-matrix(0,numvar,L)
for (i in 1:L) se.H[,i]<-sqrt(diag(psi)/(numobs*EExxy[i,i]*(1-r)))
se.psi<-diag(sqrt(2*diag(psi^2)/(numobs*(1-r))))
se.mu<- sqrt(vu/(w*numobs*(1-r)))
se.vu<- sqrt(vu^2/(w*numobs*(1-r)))
se.w<-sqrt(1/((1-w)*w*numobs*(1-r)))
std.err<-list(H=se.H,psi=se.psi,mu=se.mu,vu=se.vu,w=se.w)


output<-list(H=H,lik=likelihood,w=w,mu=mu,vu=vu,psi=psi,totni=totni,ni=ni,L=L,sigma=sigma,pqy=pqy,numvar=numvar,numobs=numobs,scaling=scaling,std.err=std.err,init=init)
invisible(output)
}
