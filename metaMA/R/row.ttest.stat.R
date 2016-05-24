`row.ttest.stat` <-
function(mat1,mat2){ 
n1<-dim(mat1)[2]
n2<-dim(mat2)[2] 
n<-n1+n2 
m1<-rowMeans(mat1,na.rm=TRUE) 
m2<-rowMeans(mat2,na.rm=TRUE) 
v1<-rowVars(mat1,na.rm=TRUE) 
v2<-rowVars(mat2,na.rm=TRUE) 
vpool<-(n1-1)/(n-2)*v1 + (n2-1)/(n-2)*v2 
tstat<-sqrt(n1*n2/n)*(m2-m1)/sqrt(vpool) 
return(tstat)}

