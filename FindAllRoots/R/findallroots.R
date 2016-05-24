##x代表未知量向量，a是系数向量,b为x的幂数,a与b是对应的
f<-function(x,a,b){
z=0
for(i in 1:length(b))
z=z+a[i]*x^(b[i])
z
}

##allroots是来求方程的所有根，a为系数,对应幂数从大到小,b为x的幂数
allroots<-function(a,b){
a1=a
b1=b
n=length(b)-1
a=a/a[1]
b=matrix(0,ncol=n,nrow=n)
for(i in 1:(n-1))
b[i,i+1]=1
for(i in 1:n)
b[n,i]=-a[n+2-i]
c=eigen(b)
print(c$values)
print("inaccuracy error")
print(f(c$values,a1,b1))
}


##x1-x2代表区间，a是系数向量,b为x的幂数,NA说明x1-x2之间没根
dichotomy<-function(x1,x2,a,b,pert = 10^(-5),n=1000,s=0.1){
x0=rep(NA,length(x1))
for(i in 1:length(x1)){
if(f(x1[i],a,b)==0)
x0[i]=x1[i]
if(f(x2[i],a,b)==0)
x0[i]=x2[i]
if(f(x1[i],a,b)!=0&f(x2[i],a,b)!=0){
x0[i]=(x1[i]+x2[i])/2
k=1
while((abs(f(x0[i],a,b))>=pert)&(k<n)){
if(f(x0[i],a,b)==0)
break
if(f(x1[i],a,b)*f(x0[i],a,b)<0)
x2[i]=x0[i]
if(f(x2[i],a,b)*f(x0[i],a,b)<0)
x1[i]=x0[i]
if(x1[i]!=x0[i]&x2[i]!=x0[i])
x2[i]=x2[i]-s
x0[i]=(x1[i]+x2[i])/2
k=k+1
if(k==1000)
x0[i]=NA
}
}
}
x0
}

