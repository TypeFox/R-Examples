pcg <-
function(A,b,M,maxiter=1e5,tol=1e-6)
{
if (missing(M)) 
{
dA<-diag(A)
dA[which(dA==0)]=1e-4
Minv=diag(1/dA,nrow=nrow(A)) 
} else Minv=solve(M)
x=rep(0,length(b))
r=b-A%*%x
z=Minv%*%(r)
p=z
iter=0
sumr2=sum(r^2)
while (sumr2>tol & iter<maxiter)
{
iter=iter+1
Ap=A%*%p
a=as.numeric((t(r)%*%z)/(t(p)%*%Ap))
x=x+a*p
r1=r-a*Ap
z1=Minv%*%r1
bet=as.numeric((t(z1)%*%r1)/(t(z)%*%r))
p=z1+bet*p
z=z1
r=r1
sumr2=sum(r^2)
}
if (iter>=maxiter) x="pcg did not converge. You may increase maxiter number."
return(x)
}
