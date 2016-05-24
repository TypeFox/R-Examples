y.v <-
function(y){
y=y[!is.na(y)]
if (length(y) == 0) return(1)
A=-a(y)*diag(length(y))
A=cbind(0,A)
A=rbind(A,1)
diag(A)=1
v=solve(A)[,length(y)+1]
v
}

