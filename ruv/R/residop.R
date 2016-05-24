residop <-
function (A, B) 
return(A - B %*% solve(t(B) %*% B) %*% t(B) %*% A)
