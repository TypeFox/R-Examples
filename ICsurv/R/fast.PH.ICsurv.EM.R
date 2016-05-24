fast.PH.ICsurv.EM <-
function(d1, d2, d3, Li, Ri, Xp, n.int, order, g0, b0, tol, t.seq, equal=TRUE){
P<-length(b0)
L<-length(g0)
N<-length(d1)


Li[d1==1]<-Ri[d1==1]
Ri[d3==1]<-Li[d3==1]
ti<-c(Li[d1==0],Ri[d3==0]) # the set of all finite endpoints of the observed intervals

if(equal==TRUE){
ti.max<-max(ti)+.00001
ti.min<-min(ti)-.00001
knots<-seq(ti.min,ti.max,length.out=(n.int+2))
}

if(equal==FALSE){
id<-seq(0,1,length.out=(n.int+2))
id<-id[-c(1,(n.int+2))]
ti.max<-max(ti)+.00001
ti.min<-min(ti)-.00001
knots<-c(ti.min,quantile(ti,id),ti.max)
}


bRi<-t(Ispline(x=Ri,order=order,knots=knots))
bLi<-t(Ispline(x=Li,order=order,knots=knots))
bt<-t(Ispline(x=t.seq,order=order,knots=knots))

###############################################
# Loading the functions needed to update the 
# parameters according to our EM algorithm

Q1<-function(b0,b1,g0,Xp,bRi,bLi,d1,d2,d3,L){
g1<-rep(-99,L)
xb0<-Xp%*%b0
xb1<-Xp%*%b1
dz<-1-exp(-(bRi%*%g0)*exp(xb0))
dw<-1-exp(-(bRi%*%g0-bLi%*%g0)*exp(xb0))
dw[d2==0]=1

EZil<-t(t(d1*bRi)*g0)*as.vector(exp(xb0)/dz)
EWil<-t(t(d2*(bRi-bLi))*g0)*as.vector(exp(xb0)/dw)
num<-EZil+(d2+d3)*EWil
den<-((d1+d2)*bRi+d3*bLi)*as.vector(exp(xb1))
g1<- apply(num,2,sum)/apply(den,2,sum)

return(g1)
}


Q2<-function(b1,b0,g0,Xp,bRi,bLi,d1,d2,d3,L){

xb0<-Xp%*%b0
xb1<-Xp%*%b1
dz<-1-exp(-(bRi%*%g0)*exp(xb0))
dw<-1-exp(-(bRi%*%g0-bLi%*%g0)*exp(xb0))
dw[d2==0]=1

EZi<-d1*(bRi%*%g0)*exp(xb0)/dz
EWi<-d2*(bRi%*%g0-bLi%*%g0)*exp(Xp%*%b0)/dw
p1<-sum((EZi+EWi)*(Xp%*%b1))

EZil<-t(t(d1*bRi)*g0)*as.vector(exp(xb0)/dz)
EWil<-t(t(d2*(bRi-bLi))*g0)*as.vector(exp(xb0)/dw)
num<-EZil+(d2+d3)*EWil
den<-((d1+d2)*bRi+d3*bLi)*as.vector(exp(xb1))
g1<- apply(num,2,sum)/apply(den,2,sum)
p2<-sum(t(EZil)*log(g1)+t(EWil)*log(g1))

p3<-sum(((d1+d2)*(bRi%*%g1) +d3*(bLi%*%g1))*exp(xb1))

res<--(p1+p2-p3)
return(res)
}

###############################################
# Lets start iterating

b1<-optim(b0, Q2, method="Nelder-Mead",b0=b0, g0=g0, Xp=Xp, bRi=bRi, bLi=bLi, d1=d1, d2=d2, d3=d3, L=L)$par
g1<-Q1(b0,b1,g0,Xp,bRi,bLi,d1,d2,d3,L)

while(max(abs(c(b0,g0)-c(b1,g1)))>tol){

b0<-b1
g0<-g1
b1<-optim(b0, Q2, method="Nelder-Mead",b0=b0, g0=g0, Xp=Xp, bRi=bRi, bLi=bLi, d1=d1, d2=d2, d3=d3, L=L)$par
g1<-Q1(b0,b1,g0,Xp,bRi,bLi,d1,d2,d3,L)
}

GRi<-bRi%*%g1
GLi<-bLi%*%g1
xb<-Xp%*%b1
ll<- sum(log((1-exp(-GRi*exp(xb)))^d1)+log((exp(-GLi*exp(xb))-exp(-GRi*exp(xb)))^d2)+log((exp(-GLi*exp(xb)))^d3))
AIC<-2*(length(b1)+length(g1))-2*ll
BIC<-(length(b1)+length(g1))*log(N)-2*ll
v<-fast.PH.Louis.ICsurv(b1,g1,bLi,bRi,d1,d2,d3,Xp)
Hessian=v

flag=is.non.singular.matrix(v)
if (flag) {var.b=solve(v)[1:P,1:P]} else {
A=v[1:P, 1:P]
B=v[1:P, (P+1):(P+L)]
C=v[(P+1):(P+L), 1:P]
D=v[(P+1):(P+L), (P+1):(P+L)]
var.b=ginv(A-B%*%ginv(D)%*%C)
} 

hz<-bt%*%g1
return(list("b"=b1, "g"=g1, "hz"=hz, "Hessian"=Hessian, "var.b"=var.b,  "flag"=flag, "ll"=ll, "AIC"=AIC, "BIC"=BIC))
}
