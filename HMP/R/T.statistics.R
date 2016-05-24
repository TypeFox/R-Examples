T.statistics <-
function(data){
T <- sum((1/colSums(data)) * apply((data-(1/sum(data))*rowSums(data) %*% t(colSums(data)))^2, 2, sum))

return(T)
}
