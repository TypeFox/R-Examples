get_z<-function(data){
s1<-sapply( 1:ncol(data) , function(y) sapply(1:nrow(data), function(x) length(unlist(data[x,y]))) )
data[s1==0]<-0
s2<-sapply( 1:ncol(data) , function(y) sapply(1:nrow(data), function(x) unlist(data[x,y])))
return(s2)}
