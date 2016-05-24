dichotomy_fun <-
function(type_kernel="n", array1,s0,s1,bandwidth,S)

  {
  maxit=1000
  tol=0.00001
  A=s0
  B=s1
  xold=s0
 for(i in 1:maxit)
  { 
  xmid=0.5*(A+B)
  ERR=abs(xmid-xold)/abs(xmid)
  if (ERR < tol) 
      xold=xmid
      else
      {
      if(functional(type_kernel, A,array1,bandwidth,S)*functional(type_kernel, xmid,array1,bandwidth,S)<0.0) B=xmid else A=xmid
      xold=xmid
      i=i+1
      }
  }
return(xold)
}
