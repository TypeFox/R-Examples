comb <-
function(m,n){
lst <- NA
if (m==0 & n==0) return(NULL)
if (m==0) return(matrix(1,1,n))
if (n==0) return(matrix(0,1,m))
k <- min(c(m,n))
for (i in 0:k){
lst <- rbind(lst, concat(comb(i,n-i), comb(m-i,i)))
}
return(lst[-1,])
}

