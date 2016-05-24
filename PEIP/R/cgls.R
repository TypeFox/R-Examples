cgls <-
function( Gmat,dee,niter)
{

 ##  
  
###  Figure out problem size.

  d1 = dim(Gmat)

  nrows = d1[1]
  ncols = d1[2]

  if(length(dee) != nrows)
    {
      print('Gmat and dee do not match in size.')
      return(NULL)
    } #

###  Setup space for the results.
  X=matrix(0,ncols,niter)
  rho=rep(0,niter)
  eta=rep(0,niter)

###  Setup for the first iteration.
  m=rep(0,ncols)
  p=rep(0,ncols)
  beta=0
  s=dee
  r= as.vector( t(Gmat)  %*% s )

###  Main loop- perform CGLS iterations.

  for (  k in 1:niter  )
    {
###  We'll precompute r'*r since it's used in several places.
      rtr= sum( r * r )

###   Update beta.
      if(k>1)
        {
          beta= rtr /  (t(prevr) %*% prevr)
        } #

###   Update p
      p= r + beta * p

###  Compute the new alpha.  To avoid doing the matrix vector
###  multiplication repeatedly, we store Gmat*p in Gp
      Gp= Gmat %*% p
      
      alpha=rtr / ( sum( Gp * Gp)) 

###  Update m.
      m=m + alpha * p

###  Update s.
      s= s - alpha *  Gp

###  Save r for the next iteration, and then update it.
      prevr=r
      r= as.vector(  t(Gmat) %*% s )

###  Store the new iterate.
      X[,k]=m
      rho[k]=Vnorm(s)
      eta[k]=Vnorm(m)
    } #



  return(list(X=X,rho=rho,eta=eta))

}
