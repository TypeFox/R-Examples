Z.statistics <-
function(data){
Num1 <- sum(data)
Num2 <- sum(apply(data, 2, function(x){
p <- sum((x-1)*x) 
return(p)
})/colSums(data)) 
Num3 <- sum(apply(data, 1, function(x){
nx <- sum(x)
b <- nx*(nx-1)
return(b)
}))
Denominator <- sqrt(2*(ncol(data)-1) * Num3)

Z <- (Num1*Num2-Num3)/Denominator

return(Z)
}
