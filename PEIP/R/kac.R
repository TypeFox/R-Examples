kac <-
function(A,b,tolx,maxiter)
{

### % First, find the size of the matrix.
  d1=dim(A);

  m = d1[1]
  n = d1[2]


### % Make a copy of A' to speed up some accesses.
  AT=t(A)

### % Setup an initial solution of all zeros.

  x = rep(0, length=n)
  iter=0;

### %  Precompute the row norms squared.
  n2=rep(0, times=m);

  
  for(i in 1:m)
    {
      n2[i]=Vnorm(AT[,i])^2
    }

### % The main loop performs iterations of Kaczmarz algorithm until 
### % maxiters is exceeded or successive iterates differ by less 
### % than tolx.  
  while(iter <= maxiter)
    {
### % Update the iteration count.
      iter=iter+1;

### %  Start the update cycle with the current solution.
      newx=x;

### %  Perform a cycle of m updates.
      for(i in 1:m)
        {
          newx=newx-(( newx  *  AT[,i]-b[i])/(n2[i])) * AT[,i]
          

### %  Check for convergence to fixed solution.
          if(  Vnorm(newx-x)/(1+Vnorm(x)) < tolx )
            {
              x=newx;
              return(x)
            }

### %  Update x for the next major iteration.
          x=newx;
        }

### % If no convergence
      
    }
   print('Max iterations exceeded.')
  return(x)
}
