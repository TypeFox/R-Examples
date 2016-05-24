PH.Louis.ICsurv <-
function(b,g,bLi,bRi,d1,d2,d3,Xp){
xb<-Xp%*%b
##############################################################
# Second partial of Q
A<-matrix(0,length(b),length(b))
B<-matrix(0,length(g),length(b))
C<-matrix(0,length(g),length(g))
# the non-zero diagonal elements of C are obsorbed in G below. 


for(i in 1:length(b)){
for(j in 1:length(b)){
A[i,j]<- -sum(((d1+d2)*bRi%*%g+d3*bLi%*%g)*exp(xb)*(Xp[,i]*Xp[,j]))
}}

for(i in 1:length(b)){
for(j in 1:length(g)){
B[j,i]<- -sum(((d1+d2)*bRi[,j]+d3*bLi[,j])*exp(xb)*Xp[,i])
}}

ci<-1-exp(-(bRi%*%g)*exp(xb))
di<-1-exp(-(bRi%*%g-bLi%*%g)*exp(xb))
di[d2==0]=1
hi<-d1*(bRi%*%g)*exp(xb)
ti<-d2*(bRi%*%g-bLi%*%g)*exp(xb)

D<-rbind(cbind(A,t(B)),cbind(B,C))

##############################################################
# Variance of the first derivative of the augmented likelihood

E<-matrix(0,length(b),length(b))
F<-matrix(0,length(g),length(b))
G<-matrix(0,length(g),length(g))

VZi<-(hi^2/ci)*(1-1/ci)+hi/ci
VWi<-(ti^2/di)*(1-1/di)+ti/di


for(i in 1:length(b)){
for(j in 1:length(b)){
E[i,j]<-sum((VZi+(d2+d3)*VWi)*Xp[,i]*Xp[,j])
}}


for(i in 1:length(b)){
for(j in 1:length(g)){
zpil<-bRi[,j]*g[j]/(bRi%*%g)
num<-(bRi%*%g-bLi%*%g)
num[d2==0]<-1
wpil<-(bRi[,j]-bLi[,j])*g[j]/num
CovZilZi<-zpil*(hi/ci)*(1+hi-hi/ci)
CovWilWi<-wpil*(ti/di)*(1+ti-ti/di)
F[j,i]<-1/(g[j])*sum((CovZilZi+(d2+d3)*CovWilWi)*Xp[,i])
}}


for(i in 1:length(g)){
for(j in 1:length(g)){

zpi<-bRi[,i]*g[i]/(bRi%*%g)
num<-(bRi%*%g-bLi%*%g)
num[d2==0]<-1
wpi<-(bRi[,i]-bLi[,i])*g[i]/num

zpj<-bRi[,j]*g[j]/(bRi%*%g)
num<-(bRi%*%g-bLi%*%g)
num[d2==0]<-1
wpj<-(bRi[,j]-bLi[,j])*g[j]/num

Covz<- zpi*zpj*(hi^2/ci)*(1-1/ci)
Covw<- wpi*wpj*(ti^2/di)*(1-1/di)

G[i,j]<- 1/(g[i]*g[j])*sum(Covz+(d2+d3)*Covw)

}}


H<-rbind(cbind(E,t(F)),cbind(F,G))

hess<--(D+H)

return(hess)
}
