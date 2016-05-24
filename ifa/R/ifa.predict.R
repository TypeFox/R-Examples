"ifa.predict" <-function(y,output,method="lms")
{
scaling<-output$scaling
if (scaling) y<-scale(y)
y<-t(y)
L<-output$L
numobs<-output$numobs
numvar<-output$numvar
ni<-output$ni
totni<-output$totni
mu<-output$mu
vu<-output$vu
sigma<-output$sigma
H<-output$H
psi<-output$psi
pqy<-output$pqy
vero<-FALSE
if (method=="bartlett") {
vero<-TRUE
x<-solve(t(H)%*%solve(psi)%*%H)%*%t(H)%*%solve(psi)%*%y}

if (method=="thompson") {
vero<-TRUE
x<-solve(diag(L)+t(H)%*%solve(psi)%*%H)%*%t(H)%*%solve(psi)%*%y }

if (method=="lms") {
vero<-TRUE
muq<-matrix(0,L,totni)
vuq<-array(0,c(L,L,totni))
pr<-totni
cont<-1
for (k in 1:L) {j<-1
                  pr<-pr/ni[k]
                  for (i in 1:totni) {muq[k,i]<-mu[k,j]
                                         vuq[k,k,i]<-vu[k,j]
                                         cont<-cont+1
                                         if (cont>pr) {if (j==ni[k]) j<-1 else (j<-j+1) 
                                                        cont<-1}
                                        }}
x<-matrix(0,L,numobs)
Aq<-array(0,c(totni,L,numvar))
bq<-matrix(0,totni,L)
sdelta<-solve(psi)

for (i in 1:totni)
    {Aq[i,,]<-sigma[,   ,i]%*%t(H)%*%sdelta
    bq[i,]<-sigma[,,i]%*%solve(vuq[,,i])%*%muq[,i]
    for (j in 1:numobs) {x[,j]<-x[,j]+(pqy[j,i]*(Aq[i,,]%*%t(t(y[,j]))+bq[i,]) )
        }}
 }
if (vero==FALSE) stop("Unrecognized method")
x<-t(x)
return(x)
}
