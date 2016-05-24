getAmoment <-
function(x,y,z=NULL){ 
# we have x (an mXn matrix), y (an n vector), and z (an n vector)
n=length(y)
mu=rowSums(x)*sum(y)/n
# line below fills in a fake z if no z has been provided
if (length(z)==0){z=rep(1,n)}
# we need to identify the classes of z
zfactor=as.factor(z)
uniquez=sort(unique(z))
K=length(unique(zfactor))
ztable=table(z)
m=dim(x)[1]
EAk=matrix(0,m,K)
EAk2=matrix(0,m,K)
EAk3=matrix(0,m,K)
EAk4=matrix(0,m,K)
# the results later will be greatly simplified if we center the y's within each
# stratum level.  This is without loss of generality, as it
# changes A by only a constant
for (k in (1:K)){
  Ik=grep(T,z==uniquez[k])
  y[Ik]=y[Ik]-mean(y[Ik])
 }

A=as.vector(x%*%y)

for (k in (1:K)){
# print(k)
 Ik=grep(T,z==uniquez[k])
 if (length(Ik)>1){
  result=getAkmoment(x[,Ik],y[Ik])
  EAk[,k]=result$final1
  EAk2[,k]=result$final2
  EAk3[,k]=result$final3
  EAk4[,k]=result$final4}
 }
# now that we have computed the per-stratum expectations,
# we need to compute the values for A, using the fact that EAk=0
EA=0
EA2=rowSums(EAk2)
EA3=rowSums(EAk3)
EA2prod=rowSums(EAk2)^2 # NOTE THAT PREVIOUS CODE DID THIS STEP STUPIDLY
EA4=rowSums(EAk4)+3*EA2prod-3*rowSums(EAk2^2)
return(list(n=n,A=A,EA=EA,EA2=EA2,EA3=EA3,EA4=EA4,mu=mu,x=x,y=y,z=z))
}



