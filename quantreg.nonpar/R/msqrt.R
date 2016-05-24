msqrt <-
function(a) 
{
  a.eig <- eigen(a);
  return(a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors));
}
