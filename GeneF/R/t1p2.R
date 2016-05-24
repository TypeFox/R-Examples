`t1p2` <-
function(v, n1, n2){
if(n1>1 && n2>1){
x <- v[1:n1]
y <- v[(n1+1):(n1+n2)]
m1 <- mean(x)
m2 <- mean(y)
sp <- (sum((x-m1)^2)+sum((y-m2)^2))/(n1+n2-2)
t <- (m2-m1)/sqrt(sp*(1/n1+1/n2))
if( is.na(t) ){
return(1)
}else{
return( 1-pt(t, n1+n2-2) )
}
}else{
return( 1 )
}
}

