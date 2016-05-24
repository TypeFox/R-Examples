noe<-function(tn,ta,tb){
n<-tn
a<-ta
b<-tb
# if (n<1) return
cgh.mat<-matrix(0,2*n+2,3)
cgh.mat<-noe.compute.cgh(n,a,b)
cv<-cgh.mat[,1]
g<-cgh.mat[,2]
h<-cgh.mat[,3]
x<-cv[-(2*n+2)]
p<-rep(0,2*n+1)
p<-noe.compute.pv(n,x)
Q0<-rep(0,2*n+2)
nQ<-rep(0,(2*n+2))
coef<-rep(0,n+1)
# Set up lower and upper bounds of recursion iterations
i1low<-h[-c(1,2)]-1
i1upp<-g[-c(2*n+1,2*n+2)]
k3low<-k2low<-k1low<-h[-c(1,2*n+2)]-1
i2low<-i1low
i2upp<-i1upp+1
# m is the primary index
Q0[1]<-1
m<-1
while (m<=2*n){
       i1<-i1low[m]
       while (i1<=i1upp[m]){
              coef[i1+1]<-1.0
              k1<-i1-1
              while (k1>=k1low[m]){
                     coef[k1+1]<-coef[k1+2]*p[m]*(k1+1)/(i1-k1)
                     k1<-k1-1}
              k2<-i1
              while (k2>=k2low[m]){
                     coef[k2+1]<-coef[k2+1]*Q0[k2+1]
                     k2<-k2-1}
              nQ[i1+1]<-0.0
              k3<-k3low[m]
              while (k3<=i1){
                     nQ[i1+1]<-nQ[i1+1]+coef[k3+1]
                     k3<-k3+1}
              i1<-i1+1}
              nQ[i2upp[m]+2]<-0.0
              i2<-i2low[m]
              while (i2<=i2upp[m]){
                     Q0[i2+1]<-nQ[i2+1]
                     i2<-i2+1}
       m<-m+1}
ans<-nQ[n+1]} 
