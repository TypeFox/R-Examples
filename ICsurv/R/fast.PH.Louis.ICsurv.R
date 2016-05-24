fast.PH.Louis.ICsurv <-
function(b,g,bLi,bRi,d1,d2,d3,Xp){

xb<-Xp%*%b

A<--t(as.vector(((d1+d2)*bRi%*%g+d3*bLi%*%g)*exp(xb))*Xp)%*%Xp
B<--t(((d1+d2)*bRi+d3*bLi)*as.vector(exp(xb)))%*%Xp
C<-matrix(0,length(g),length(g))
D<-rbind(cbind(A,t(B)),cbind(B,C))

ci<-1-exp(-(bRi%*%g)*exp(xb))
di<-1-exp(-(bRi%*%g-bLi%*%g)*exp(xb))
di[d2==0]=1
hi<-d1*(bRi%*%g)*exp(xb)
ti<-d2*(bRi%*%g-bLi%*%g)*exp(xb)
VZi<-(hi^2/ci)*(1-1/ci)+hi/ci
VWi<-(ti^2/di)*(1-1/di)+ti/di

E<-t(as.vector((VZi+(d2+d3)*VWi))*Xp)%*%Xp

zpil<-t(t(bRi)*g)/as.vector(bRi%*%g)
num<-(bRi%*%g-bLi%*%g)
num[d2==0]<-1
wpil<-t(t(bRi-bLi)*g)/as.vector(num)
CovZilZi<-zpil*as.vector((hi/ci)*(1+hi-hi/ci))
CovWilWi<-wpil*as.vector((ti/di)*(1+ti-ti/di))
F<-t(CovZilZi+(d2+d3)*CovWilWi)%*%Xp/g

Covz<-t(zpil*as.vector((hi^2/ci)*(1-1/ci)))%*%zpil
Covw<-t(wpil*as.vector((d2+d3)*(ti^2/di)*(1-1/di)))%*%wpil
G<-t((Covz+Covw)/g)/g

H<-rbind(cbind(E,t(F)),cbind(F,G))

hess<--(D+H)
return(hess)
}
