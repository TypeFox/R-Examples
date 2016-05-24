change<-function(levset){
#
#
len<-length(levset)
m<-sum(levset)
rindeksit<-matrix(0,m,1)
j<-1
for (i in 1:len){
    if (levset[i]==1){
       rindeksit[j]<-i
       j<-j+1
    }
}
return(rindeksit)
}

