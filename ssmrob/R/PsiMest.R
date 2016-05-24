PsiMest <-
function(x, y, beta, sigma, t.c = 1.345, weight = 1)  
{ 
  if((y-t(as.vector(x))%*%as.vector(beta))/sigma<(-t.c))  
    return(-t.c*x*weight)
  else if((y-t(as.vector(x))%*%as.vector(beta))/sigma>t.c) return(t.c*x*weight)
  else return((y-t(as.vector(x))%*%as.vector(beta))/sigma*x*weight)
}
