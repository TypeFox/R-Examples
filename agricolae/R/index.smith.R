`index.smith` <-
function(data,...) {
A<-as.matrix(data)
n<-nrow(A)
m<-ncol(A)
TC <- sum(A)^2/(n*m)
media <- mean(A)
m1<-1:m
n1<-1:n
k0<-m1[m/m1==m%/%m1] # para las columnas
r0<-n1[n/n1==n%/%n1] # para las filas
lk0 <-length(k0)
lr0 <-length(r0)
k0 <- k0[-lk0]
r0 <- r0[-lr0]
lk0 <-length(k0)
lr0 <-length(r0)
x<-rep(0,lk0*lr0);V<-x; CV<-x;Width<-x;Length<-x;l<-0
for (k in k0) {
for (p in r0) {
l<-l+1
ss<-0
k3<-0;k4<-0
for (i in 1:(n/p)) {
k3<- k4+1
k4<- k3+p-1
k1<-0;k2<-0
for (j in 1: (m/k)) {
k1<- k2+1
k2 <-k1+k-1
ss<-ss+(sum(A[k3:k4,k1:k2]))^2/(k*p)
}
}

V[l]<-ss-TC
#V[l]<-V[l]/(k*p)^2
Length[l]<-p
Width[l]<-k

x[l] <- k*p
V[l]<- V[l]/(n*m-1)
CV[l]<-sqrt(V[l])*100/media
}
}
z<-log(x[-1])
y<-log(V[1]/V[-1])
model <-lm(y ~ 0+z)
b <- as.numeric(coef(model))
tabla <-cbind(Size=x,Width,Length,plots=n*m/x,Vx=V,CV=round(CV,1))
uniformity<-tabla[order(tabla[,1]),]
model <- lm(CV ~ I(log(x)))
coeff <- coef(model)
size<-1:max(x)
cv<- coeff[1]+coeff[2]*log(size)
plot(size,cv,...)
points(x,CV)
return(list(model=model,uniformity=uniformity))
}

