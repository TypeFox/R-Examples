MmatrM <-
function(x, y, beta, sigma, t.c = 1.345, weights = 1)  
{
  (t(x)%*%diag(drop(weights*ifelse(abs(y-(x)%*%as.vector(beta))/sigma<t.c,1,0)),length(x[,1]),length(x[,1]))%*%x)/sigma 
}
