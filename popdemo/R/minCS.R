minCS <-
function(A){
if(length(dim(A))!=2) stop("A must be a matrix")
order<-dim(A)[2]
add<-(1:order)
for(i in 1:order){
    add[i]<-sum(A[,i])
}
return(min(add))
}

