"balancedcluster" <-
function(X,m,cluster,selection=1,comment=TRUE,method=1)
{
cluster=cleanstrata(cluster)
if(comment==TRUE) cat("\nSELECTION OF A SAMPLE OF CLUSTERS\n")
p=dim(X)[2]
N=dim(X)[1]
H=max(cluster)
XC=array(0,c(H,p))
Ni=rep(0,times=H)
for(h in 1:H)
   { Ni[h]=sum(as.integer(cluster==h)) 
    for(j in 1:p)  XC[h,j]=sum(X[cluster==h,j]) 
   }
if(selection==1) pik=inclusionprobabilities(Ni,m) else pik=rep(m/H,times=H)
s=samplecube(cbind(pik,XC),pik,1,comment,method) 
res=array(0,c(N,2))
for(h in 1:H) 
    {
     res[cluster==h,1]=s[h] 
     res[cluster==h,2]=pik[h]
     }
res
}


