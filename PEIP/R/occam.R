occam <-
function( afun,ajac,L,d,m0,delta)
{
  #############   for the condition number
####  library(pracma)
####  library(Matrix)

  m=m0
  oldm=rep(0,length(m))
  iter=0
#############    need this to start
  mchi2 = delta^2*1.01
       func =  get(afun)
      jacob =  get(ajac)

###  while we have not converged sufficiently or the data misfit is higher than 
###  allowed keep iterating
  while((Vnorm(oldm-m)/Vnorm(m) >5.0e-3) | (mchi2 > delta^2*1.01) )
    {
###  only allow 30 iterations
      iter=iter+1
      cat(paste("Iteration:", iter), sep="\n")
      
      if(iter > 30)
        {
          return(m)
        }

###  store the old mode to test for convergance
      oldm=m

###  get the current data that would be generated and the jacobian

      G=func(m)
      J=jacob(m)

###  get the dhat that is in equation 10.14
      dhat=d-G+  J %*% m 
      
###  This is a simple brute force way to do the line search.  Much more
###  sophisticated methods are available.  Note: we've restricted the line
###  search to the range from 1.0e-20 to 1.  This seems to work well in
###  practice, but might need to be adjusted for a particular problem.

      alphas=RSEIS::logspace(-20,0,100)

      chis = rep(0, length=100)
      
      
      for(i in 1:100 )
        {
          M=t(J) %*% J +alphas[i]^2 * t(L) %*% L

###  if M is not terribly conditioned
          if( Matrix::condest(M)$est < 1.0e15  )
            {
              m= as.vector( solve( t(J)%*%J+alphas[i]^2*t(L)%*%L  ) %*%   t(J) %*% dhat )


              
###  m=Ainv( t(J)%*%J+alphas[i]^2*t(L)%*%L  , as.vector( t(J) %*% dhat )    )   

              
###  store the associated data misfit
              chis[i]=Vnorm(func(m)-d)^2
            }
          else
            {
###  M behaves poorly enough it should not be used
              chis[i]=+Inf
            }
        }

      Y=min(chis)
      I = which.min(chis)
      

      if (Y > delta^2)
        {
### disp('Improving Chi^2')
          alpha=alphas[I[1]]
        }
      else
        {
### disp('Smoothing m')
          I=which( chis<=delta^2 )
          alpha=alphas[max(I)]
          
        }

###  store the new model and misfit
      m=Ainv( (t(J) %*% J +alpha^2*t(L)%*%L) ,  as.vector( t(J) %*% dhat)  )
      
      mchi2=Vnorm( func(m)-d )^2
      
    }


  return(m)

}
