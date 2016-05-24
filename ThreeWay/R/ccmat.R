ccmat <-
function(A,B){

s=dim(A)
if (s[1]<nrow(B)){
	s[1]=nrow(B)
}
if (s[2]<ncol(B)){
	s[2]=ncol(B)
}
mat=matrix(0,s[1],s[2]*2)
j=1
k=1
for (i in 1:ncol(mat)){
	if ((i%%2)==1){
		mat[1:s[1],i]=A[1:s[1],j]
		j=j+1
	} else{
		mat[1:s[1],i]=B[1:s[1],k]
		k=k+1
	}
}
return(mat)
}
