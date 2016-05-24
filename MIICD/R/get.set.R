
get.set <- function( k , data){
data <- data[,c('left','right')]
df1<-data.frame( data , r = runif(nrow(data)))
w <- with( df1 , left +   (right - left)  * r ) 
if1<-ifelse(w == Inf , df1$left , w)
return(if1)
}


