is.even <-
function(x){
start <- 1
end <- length(x)+1
while(start<end){
y <- x[start]

test1 <- y/2
test2 <- floor(test1)
if(test1 == test2) {
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

