`get.numbers` <-
function(x){
n <- tapply(x, x, length)
l <- double(length(n)+1)
for(i in 1:length(n)){
l[i+1] <- l[i]+n[i]
}
return(l)
}

