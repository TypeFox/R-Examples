disS=function(p){
# %%%Compute density based Silhouette information (DBS)
# %%%%%%INPUT
# %p probability matrix
# %%%%%%OUTPUT
# %ds distances for silhouette
  
n=nrow(p)
nc=ncol(p)
m=matrix(0,n,1)
for( i in 1:n){
  m[i]=max(p[i,])}
pm=max.col(p)
pt=matrix(0,n,(nc-1))
for( i in 1:n){
  c=0
  for(j in 1:nc){
    if(j!=pm[i]){
      c=c+1
      pt[i,c]=p[i,j]
    }
  }  
}
nu=matrix(0,n,1)
for( i in 1:n){
nu[i]=log(m[i]/max(pt[i,]))}
de2=max(nu)
ds=nu/de2
return(ds)}