## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

subspace_min <- function(x, l, u, x_cp, d, G){
n=length(x)
Z=diag(n)
fixed=(x_cp<=l+1e-8 | x_cp>=u-1e8)
if(all(fixed)){
    x=x_cp
    return(x)
}
Z=Z[,!fixed]
rgc=t(Z)%*%(d+G%*%(x_cp-x))
rB=t(Z)%*%G%*%Z
d[!fixed]=solve(rB,rgc)
d[!fixed]=-d[!fixed]
alpha=1
temp1=alpha
for(i in find_matlab(!fixed)){
    dk=d[i]
    if(dk<0){
        temp2=l[i]-x_cp[i]
        if(temp2>=0){
            temp1=0
        } else {
            if(dk*alpha < temp2) temp1=temp2/dk
        } 
    } else {
        temp2=u[i]-x_cp[i]
        if(temp1<=0){
            temp1=0
        } else {
            if(dk*alpha>temp2) temp1=temp2/dk
        }
    }
    alpha = min(temp1,alpha)
}
x=x_cp+Z*alpha*d[!fixed]
return( x)
}
