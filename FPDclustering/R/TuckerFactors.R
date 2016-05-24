TuckerFactors<- function(data=NULL,  nc=2){

# %heuristic approach to cope with the choice of Tucker 3 factors.
# %we choose the minimum number of factors that corresponds to a significant value 
# %of the explained variability.
# %
# %%%%%%%%%INPUT%%%%%%%%%%
# %dataO dataset
# % nu number of clusters
# %mi minimum number of factors
# %%%%%%%%%ouPUT%%%%%%%%%%
# %expleined expleined variability
# 
# 
# 
# %mi=1;

n=nrow(data)
J=ncol(data)


#%Computation of distance array G
pd=PDclust(data,nc)
c=pd$centers

ddt=matrix(0,n,nc*J)
dd=matrix(0,n,J)
for( j in 1:nc){
  for( i in 1:n){
    v=t(as.vector(abs(data[i,]-c[j,])))
    dd[i,]=v
  }
  ddt[,((1+J*(j-1)):(J*(j-1)+J))]=dd
}

m=10
if(J<10){ m=J}

nfc=nc
ini=2#the minimum number of factor for variables is 2
if( nc>2){
nfc=nc-1# if there are more than two clusters the number of factor for clusters is nc-1;
}
expleined=matrix(0,(m-(ini-1)),19)
  niter=(m-ini+1)*19
  cont=0
for( nf in ini:m){
  cont=cont+1
  cat(paste(round(cont/niter*100)),"% of the process completed", fill = TRUE)
explold=0
for( nu in 2:20){   
  cont=cont+1
 
tuk3=T3funcrep(ddt, n, J, nc, nu, nf,nfc, start=0,conv=0.1)#tucker3decomposition
expl=tuk3$fp
expleined[(nf-(ini-1)),(nu-1)]=expl#%explaned variability is stored 

if( expl<explold){
expleined[(nf-(ini-1)),(nu-1)]=explold}
expl=explold
}
explold=expl
}
len=nrow(expleined)
leg=c('10 factors','9 factors','8 factors','7 factors','6 factors','5 factors','4 factors','3 factors','2 factors')
leg=leg[(9-len+1):9]
plot((2:20),expleined[len,],type='l',main='Explained variability',ylab='Explained variability', xlab='number of factors for units',col=len)
  axis(1,at=2:20)
  for(i in (len-1):1){
  matplot((2:20),expleined[i,],type='l',add=T,col=(i))
}
legend(15,max(expleined),leg,text.col=(len:1))#(expleined[10,3])
return(expleined)
}