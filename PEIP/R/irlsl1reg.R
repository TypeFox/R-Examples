irlsl1reg <-
function( G,d,L,alpha,maxiter=100,tolx=1.0e-4,tolr=1.0e-6)
{

  ##

 
###  unchanging constants in the system that is solved repeatedly
  GTG=t(G) %*% G
  GTd=t(G) %*% d

###  Start with an initial unweighted solution.
  m= Ainv(  (2*GTG+alpha*(t(L) %*%L)) , (2*GTd) )

  iter=1

###  iterate until maxiter or we converge and return
  while(iter < maxiter)
    {
      iter=iter+1

###   get the magnitude of Lm, but don't let any element be less than tolr
      absLm=abs(L %*% m)
      absLm[absLm<tolr]=tolr

###  build the diagonal weighting matrix for this iteration
      R=diag(1/as.vector(absLm) )

      mold=m

      KADD = alpha * t(L) %*% R %*% L
   ###   print(dim(KADD))
      
###  get the new iterate and check for convergance
      m = Ainv( (2*GTG + KADD  ) , (2*GTd) )


      if(Vnorm(m-mold)/(1+Vnorm(mold)) < tolx)
        {
          mreg=m
          return(mreg)
        } 
    } 

###  Give a warning, if desired, but return best solution.
### warning('irlslreg1 maximum iterations exceeded.')
  mreg=m

  return(mreg)

}
