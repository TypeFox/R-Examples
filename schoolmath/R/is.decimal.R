is.decimal <-
function(x){
start <- 1
end <- length(x)+1
while(start<end){
y <- x[start]

test <- floor(y)
if(y==test){
if(start==1){
result=FALSE 
}else{
result<- c(result,FALSE)
}

}else{
if(start==1){
result=TRUE
}else{
result <- c(result,TRUE)
}
}
start <- start+1
}

return(result)
}

