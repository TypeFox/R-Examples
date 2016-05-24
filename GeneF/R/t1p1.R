`t1p1` <-
function(v, n){
if(n>1){
m <- mean(v)
s <- var(v)
t <- m/sqrt(s/n)
if( is.na(t) ){
return(1)
}else{
return( 1-pt(t, n-1) )
}
}else{
return(1)
}
}

