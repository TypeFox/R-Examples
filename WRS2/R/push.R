push<-function(mat){
#
# For every column of mat, move entry down 1
#
matn<-matrix(NA,nrow=nrow(mat),ncol=ncol(mat))
Jm<-nrow(mat)-1
for (k in 1:ncol(mat)){
temp<-mat[,k]
vec<-0
vec[2:nrow(mat)]<-temp[1:Jm]
matn[,k]<-vec
}
matn
}
