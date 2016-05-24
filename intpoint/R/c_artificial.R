c_artificial <-
function(c,A){
col<-ncol(A)
l<-length(c)
cn<-array(0,c(col,1))
for(i in 1:l)
cn[i]<-c[i]
cn[col]<-1000000
return(cn)}
