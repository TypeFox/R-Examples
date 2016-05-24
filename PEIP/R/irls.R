irls <-
function(A,b,tolr,tolx,p,maxiter)
{
  Mnorm <-function(X)
    {
### the norm function in matlab
      ##  returns the largest singular value of the matrix or vector
      ##  here we create a function to duplicate this bizare behavior    
      s = svd(X, 0, 0)

      return(max(s$d))

    }


   eps = .Machine$eps

###% Find the size of the matrix A.
###[m,n]=size(A);

  d = dim(A)

  m = d[1]
  n = d[2]

###% Start the first iteration with R=I, and x=A\b (the least
###% squares solution)
  R=diag(rep(1,times=m) );

  x = as.vector( solve(t(A) %*% A, tol=eps) %*% t(A) %*% b )

### x=A\b;

###% Now loop up to maxiter iterations
  iter=1;
  while(iter <= maxiter)
    {
      iter=iter+1;

###% compute the current residual
      r= as.vector( A%*%x - b )

###% for each row adjust the weighting factor R based on the residual
      for(i in 1:m)
        {
          if (abs(r[i]) < tolr)
            {  r[i]=abs(tolr)^(p-2); }
          else
            {
              r[i]=abs(r[i])^(p-2);
            }
        }

###% put the weighting factors into R
      R=diag(r);

###% find the solution to the weighted problem

      newA = (t(A) %*% R %*% A)
      newb =t(A) %*% R %*% b
      newx=as.vector( solve(t(newA) %*% newA, tol=eps )%*% t(newA) %*% newb )

###%check for convergence
  if(Mnorm(newx-x)/(1+Mnorm(x)) < tolx)
    {
      x = newx;
      return;
    }
  else
    {
      x=newx;
    }
}

return(x)

}
