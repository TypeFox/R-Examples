Dx <-
function(x){
l<-length(x)
dx<-array(0,c(l,l))
for (i in 1:l)
dx[i,i]<-x[i]
return(dx)}
