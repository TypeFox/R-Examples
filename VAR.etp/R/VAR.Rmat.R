VAR.Rmat <-
function(p,k,restrict,type="const")
{
info = restrict[,1:3,drop=F]

if(type=="none") add <- 0
if(type=="const") add <- 1
if(type=="const+trend") add <-2

M = p*(k^2)+add*k - nrow(info)
Rmat = diag(M); 

rmat = matrix(0,nrow=p*(k^2)+add*k)

mat1 <- rep(1:k,k)
mat2 <- rep(1:k,each=k)
mat3 <- rep(1:p,each=k^2)
position <- cbind(mat3,mat1,mat2)

tem1=numeric()
for (i in 1:nrow(position)){
for (j in 1:nrow(info)){
tem = prod(as.numeric(position[i,] == info[j,]))
if (tem == 1) tem1=c(tem1,i)}}

for(i in 1:length(tem1)){
index = tem1[i]
if(index == 1) Rmat <- rbind(matrix(0,ncol=M),Rmat)
else{
if(index <= nrow(Rmat)) Rmat <- rbind(Rmat[1:(index-1),,drop=F],matrix(0,ncol=M),Rmat[index:nrow(Rmat),,drop=F])
if(index > nrow(Rmat)) Rmat <- rbind(Rmat,matrix(0,ncol=M)) }
}
if(ncol(restrict) == 4) rmat[tem1,] = restrict[,4]
index1 = rowSums(Rmat)
index2 = 1
Cmat=matrix(0,nrow=nrow(restrict),ncol=p*k^2+add*k)
for (i in 1:length(index1))
if (index1[i] == 0) {Cmat[index2,i] = 1; index2=index2+1}
cmat = matrix(0,nrow=nrow(restrict))
if(ncol(restrict) == 4) cmat = matrix(restrict[,4])

return(list(Rmat=Rmat,rvec=rmat,Cmat=Cmat,cvec=cmat))
}
