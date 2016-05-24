is.real.positive <-
function(x){
start <- 1
end <- length(x)+1
while(start<end){
y <- x[start]

if(y>=0) {
if(start==1){
result=TRUE 
}else{
result<- c(result,TRUE)
}

}else{
if(start==1){
result=FALSE
}else{
result <- c(result,FALSE)
}
}
start <- start+1
}
return(result)
}

