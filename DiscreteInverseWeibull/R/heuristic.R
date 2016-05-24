heuristic <-
function(x, beta1=1, z=0.1, r=0.1, Leps=0.01)
{
vbeta<-numeric(3)
vq<-numeric(3)
L<-numeric(3)

beta1<-1

betastar0<-beta1
betastar<-beta1+1
Lold<-Inf
Lnew<-0

vbeta[2]<-betastar0
vbeta[1]<-vbeta[2]-z
vbeta[3]<-vbeta[2]+z

res1<-optimize(f = loglikediw, x = x, beta=vbeta[1], interval=c(0,1))
res2<-optimize(f = loglikediw, x = x, beta=vbeta[2], interval=c(0,1))
res3<-optimize(f = loglikediw, x = x, beta=vbeta[3], interval=c(0,1))

vq[1]<-res1$minimum
vq[2]<-res2$minimum
vq[3]<-res3$minimum

L[1]<-res1$objective
L[2]<-res2$objective
L[3]<-res3$objective

betastar<-vbeta[which.min(L)]

while(betastar!=betastar0 | abs(Lnew-Lold)>Leps)
{
Lold<-Lnew
betastar0<-betastar

vbeta[2]<-betastar
vbeta[1]<-vbeta[2]-z
vbeta[3]<-vbeta[2]+z

res1<-optimize(f = loglikediw, x = x, beta=vbeta[1], interval=c(0,1))
res2<-optimize(f = loglikediw, x = x, beta=vbeta[2], interval=c(0,1))
res3<-optimize(f = loglikediw, x = x, beta=vbeta[3], interval=c(0,1))

vq[2]<-res2$minimum

L[1]<-res1$objective
L[2]<-res2$objective
L[3]<-res3$objective

betastar<-vbeta[which.min(L)]

if(betastar==betastar0)
{
z<-z*r
Lnew<-min(L)
}
}

c(vq[2],vbeta[2])

}

