quadlin <-
function(Q,A,b)
  {

### % First, find lambda.
    lambda=solve( (A %*% solve(Q) %*% t(A) ), b)

### % Now, x.
    x=solve(Q)%*%t(A) %*% lambda;

    return(list(x=x,lambda=lambda))

  }
