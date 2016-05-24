## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

cauchy_point <- function( x, l, u, g, B ){
#finds the cauchy point for q(x)=x'Gx+x'd s$t. l<=x<=u

#g=G*x+d #gradient of q(x)
n=length(x)
t=zeros(n,1)
d=zeros(n,1)
for(i in 1:n){
    if(g[i]<0){
        t[i]=(x[i]-u[i])/g[i]
    }
    if(g[i]>0){
        t[i]=(x[i]-l[i])/g[i]
    }
    if(g[i]==0){
        t[i]=Inf
    }
    if(t[i]==0){
        d[i]=0
    } else {
        d[i]=-g[i]
    }
}
ts=unique(t)
ts=ts[ts!=0]
ts=sort(ts)

df=t(g)%*%d
d2f=t(d)%*%B%*%d
dt_min=-df/d2f
t_old=0
i=1
z=zeros(n,1)
while(i<=length(ts) && dt_min>=(ts[i]-t_old)){
     ind=ts[i]<t
     d[!ind]=0
     z=z+(ts[i]-t_old)*d
     df=t(g)%*%d+t(d)%*%B%*%z
     d2f=t(d)%*%B%*%d
     dt_min=-df/d2f
     t_old=ts[i]
     i=i+1
}
dt_min=max(dt_min,0,na.rm=TRUE)
t_old=t_old+dt_min
x_cp=x-t_old*g
temp=x-t*g
x_cp[t_old>t]=temp[t_old>t]

#x_cp=x+t_old*d
return(  x_cp  )
}
