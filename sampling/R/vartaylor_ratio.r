vartaylor_ratio=function(Ys,Xs,pikls)
{
if (any(is.na(pikls))) 
        stop("there are missing values in pikls")
if(nrow(pikls)!=ncol(pikls))
        stop("pikls is not a square matrix")
if (any(is.na(Ys))) 
        stop("there are missing values in y")
if (any(is.na(Xs))) 
        stop("there are missing values in x")
if (length(Ys) != nrow(pikls) | length(Xs) != nrow(pikls) | length(Xs) != length(Ys)) 
      stop("y, x and pikls have different sizes")
pik=diag(pikls)
n=length(pik)
xhat=sum(Xs/pik)
yhat=sum(Ys/pik)
r=yhat/xhat
z=(Ys-r*Xs)/xhat
delta=matrix(0,nrow=nrow(pikls),ncol=ncol(pikls))
for(i in 1:(n-1))
 {for(j in (i+1):n)
      delta[i,j]=delta[j,i]=(1-pik[i]*pik[j]/pikls[i,j])*z[i]*z[j]/(pik[i]*pik[j])
  delta[i,i]=(1-pik[i])*z[i]^2/pik[i]^2}
delta[n,n]=(1-pik[n])*z[n]^2/pik[n]^2
list(ratio=r, estvar=sum(delta))
}



