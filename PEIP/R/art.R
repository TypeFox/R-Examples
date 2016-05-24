art <-
function(A,b,tolx,maxiter)
  {

### % Alpha is a damping factor.  If alpha<1, then we won't take full steps
### % in the ART direction.  Using a smaller value of alpha (say alpha=.75)
### % can help with convergence on some problems.  
    alpha=1.0;

### % First, get the size of A.

    d1=dim(A);

    m = d1[1]
    n = d1[2]

### % In the A1 array, we convert all nonzero entries in A to 1.
    A1 = A
    A1[A>0] = 1

### % Get transposed copies of the arrays for faster access.
    AP=t(A)

    A1P=t(A1)

### %  Precompute N(i) and L(i) factors.
    N=rep(0, m);
    L=rep(0, m);
    for(i in 1:m)
      {
        N[i]=sum(A1[i,]);
        L[i]=sum(A[i, ]);
      }

### % Start with the zero solution.
    x=rep(0, n);

### % Start the iteration count at 0.
    iter=0;

### % Now, the main loop.
    while(TRUE)
      {
### % Check to make sure that we haven't exceeded maxiter.
        iter=iter+1;
        if(iter > maxiter)
          {
            print('Max iterations exceeded.');
            x=newx;
            return(x)
          }
        

### % Start the next round of updates with the current solution.
        newx=x;
        
### % Now, update each of the m constraints.
        for(i in 1:m)
          {
### %  Compute the weighted sum for constraint i.
            q=A1P[,i] * newx;

### % We use the more accurate formula for delta.
            delta=b[i]/L[i]-q/N[i];

### % This alternative formula is less accurate and doesn't work nearly as well.
### %
### %    delta=(b(i)-q)/N(i);
### %

### % Now do the update.
            newx=newx+alpha*delta*A1P[,i]
          }

### % Check for convergence
        if(Vnorm(newx-x)/(1+Vnorm(x)) < tolx)
          {
            x=newx;
            return(x)
          }

### % No convergence, so setup for the next ART iteration.
        x=newx;
        return(x)

      }

     return(x)
  }
