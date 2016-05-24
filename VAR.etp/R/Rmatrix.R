Rmatrix <-
function(p,k,type=1,index=0){

if(type==1) {Rmat=diag(k*p); rvec=matrix(0,nrow=k*p)}
if(type==2) {a=matrix(1,ncol=p); b=diag(k); Rmat= b %x% a; rvec=matrix(0,nrow=k)}

if (type==1 & index ==0) return(list(Rmat=Rmat,rvec=rvec))
if (type==1 & index > 0) 
{
indexmat=matrix(0,nrow=k,ncol=p)
index1=1:p
for (i in 1:k){
indexmat[i,]=index1
index1=index1+p}
return(list(Rmat=Rmat[indexmat[index,],,drop=F], rvec=rvec[indexmat[index,],,drop=F]))}

if (type==2 & index ==0) return(list(Rmat=Rmat,rvec=rvec))
if (type==2 & index > 0) return(list(Rmat=Rmat[index,,drop="F"], rvec=rvec[index]));
}
