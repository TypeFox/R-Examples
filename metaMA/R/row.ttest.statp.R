`row.ttest.statp` <-
function(mat){ 
m<-rowMeans(mat,na.rm=TRUE) 
sd<-sqrt(rowVars(mat,na.rm=TRUE))  
tstat<-m/(sd*sqrt(1/dim(mat)[2])) 
return(tstat)}

