logratio <- function(n,sc,lambda1, lambda2){
Q <- (1-exp(-lambda1))/(1 - exp(-lambda2))
if ( Q < 1)
{lr <-  log(n-1) + sc*log(Q) + log(1 - Q) - log(1 - Q^n)}
else
{lr <-  log(n-1) + sc*log(Q) + log(Q - 1) - log(Q^n - 1)}
} 
