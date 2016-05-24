titv <-
function(DNAbin){
mat<-as.matrix(DNAbin)
res<-matrix(NA, ncol=dim(mat)[1], nrow=dim(mat)[1], dimnames=list(x=names(DNAbin), y=names(DNAbin)))
for(i in 1:(dim(mat)[1] - 1)){
for(j in (i+1):dim(mat)[1]){
vec<-as.numeric(mat[i,])+as.numeric(mat[j,])-8
res[j,i]<-sum(!is.na(match(vec,c(200,56))))#Transitions
res[i,j]<-sum(!is.na(match(vec,c(152,168,88,104))))#Transversions
}
}
res
}

