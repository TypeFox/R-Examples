LICurvature<-function(ini = NA,X,Xstar,y,n,p, ...)UseMethod("LICurvature")
class(LICurvature)<-"LICurvature"
LICurvature.default <- function(ini = NA, X,Xstar, y,n,p,...)
{
f <- function(ini = NA, X,Xstar, y,n,p) {
p=length(X[1,])
n = nrow(X)
my = matrix(0,n,p)
l=vector("numeric",n)
l=as.vector(l)
y <- as.vector(y)
ini <- as.vector(ini)
Xstar=cbind(1,X)
Xstar <- as.matrix(Xstar)
muy=numeric(n)
for(i in 1:n){
for(j in 1:p){
my[i,j]=ini[0:p][[j]]*Xstar[i, ][[j]]
}}
for (i in 1:n) {
muy[i] <- sum(my[i,])

l[i] <-  log(dnorm(y[i], muy[i],ini[p+1]))
}

Like=sum(l)
;-Like
}
ini=c(0,rep(1,p))
ml = nlminb(ini, f, X = X,Xstar=Xstar, y=y,n=n,p=p, lower = c(rep(-Inf,
p), 0), upper = c(rep(Inf,p+1)), hessian = T)
phat=as.matrix(ml$par[0:p+1])
yhat=Xstar%*%phat
res=y-yhat
res=as.vector(res)
E=diag(res)
Px=Xstar%*%ginv((t(Xstar)%*%Xstar))%*%t(Xstar)
esq=(res)^2
esq=matrix(esq,n,1)
z1=t(Xstar)%*%E/ml$par[p+1]
z2=t(esq)/2*(ml$par[p+1])^2
Delta=rbind(z1,z2)
v1=t(Xstar)%*%Xstar/ml$par[p+1]
v2=n/2*(ml$par[p+1])^2
z=matrix(rep(0),(p+1))
m1=rbind(v1,t(z))
m2=rbind(z,v2)
Lzegond=(-1)*cbind(m1,m2)
Fzegond=t(Delta)%*%ginv(Lzegond)%*%(Delta)
eig=eigen(Fzegond)
eigenval=eig$values
maxval=max(eigenval)
eigvec=eig$vectors
maxvec=eigvec[,which((eigenval)==maxval)]
Q=t(maxvec)%*%(E%*%Px%*%E+esq%*%t(esq)/(2*n*ml$par[p+1]))%*%maxvec
C=2*abs(Q)/(ml$par[p+1])

r <- list(call = ml, lmax=maxvec,Cmax=C)
 r$call <- match.call()
    class(r) <- "LICurvature"
    r

}

#package.skeleton(name="LICurvature", code_files="I:\\LICurvature.R")
