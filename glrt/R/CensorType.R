CensorType <-
function(A, inf = Inf) 
{
cens = as.integer(rep(2, nrow(A)))  
if(any(A[,1] == 0 && A[,2] == inf))
{
warning("Observation: (0, inf)!")
return(NULL)
}
else
{
cens[which(A[,1] == 0)] = 1  
cens[which(A[,2] == inf)] = 3  

cens[which(A[,1] == A[,2])] = 4 
}
cens
}

