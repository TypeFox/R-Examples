# This is for an arbitrary matrix D. It is really just a wrapper
# for the hard work that is done by the functions dualpathWide,
# dualpathWideSparse, or dualpathTall.

dualpath <- function(y, D, approx=FALSE, maxsteps=2000, minlam=0,
                     rtol=1e-7, btol=1e-7, verbose=FALSE) {
  m = nrow(D)
  n = ncol(D)
  
  # If D has no rows, there is nothing to do
  if (m==0) {
    return(list(lambda=Inf,beta=as.matrix(y),fit=as.matrix(y),
                u=as.matrix(rep(0,m)),hit=TRUE,df=n,y=y,
                completepath=TRUE,bls=y))
  }
  
  # Get rid of entirely zero columns in D
  n0 = n
  y0 = y
  D0 = D
  j = which(colSums(D!=0)>0)
  y = y[j]
  D = D[,j,drop=FALSE]
  n = length(j)
  coldif = n0-n
  
  # If there are no columns left, there is nothing to do
  if (n==0) {
    return(list(lambda=Inf,beta=as.matrix(y0),fit=as.matrix(y0),
                u=as.matrix(rep(0,m)),hit=TRUE,df=n0,y=y0,
                completepath=TRUE,bls=y))
  }

  # This is a flag for whether we should be treating D as sparse
  sparse = c(attributes(class(D))$package,"")[[1]] == "Matrix"

  # Wide case
  if (m <= n) {

    # Compute a QR factorization of D^T
    x = qr(t(D))                
    R = qr.R(x) # m x m

    # Check if rank(D)=m, which we need for the sparse wide algorithm
    if (all(abs(diag(R)) >= rtol)) {
      # Sparse algorithm
      if (sparse) {
        out = dualpathWideSparse(y,D,x,approx,maxsteps,minlam,rtol,btol,verbose)
      }
      # Dense algorithm
      else {
        Q1 = qr.Q(x)               # n x m 
        Q2 = matrix(nrow=n,ncol=0) # n x 0
        out = dualpathWide(y,D,Q1,Q2,R,approx,maxsteps,minlam,rtol,btol,verbose)
      }
      
      # Construct beta, fit, y, bls, while accounting for the fact that
      # we may have had zero columns in D 
      out$df = out$df + coldif
      beta = matrix(y0,n0,length(out$lambda))
      beta[j,] = as.matrix(y0[j] - t(D0[,j])%*%out$u)
      colnames(beta) = colnames(out$u)
      out$beta = beta
      out$fit = beta
      out$y = y0
      out$bls = y0

      # Add to pathobjs component
      out$pathobjs$n0 = n0
      out$pathobjs$y0 = y0
      out$pathobjs$j = j
      out$pathobjs$D0 = D0
      out$pathobjs$coldif = coldif
      return(out)
    }

    else if (sparse) {
      warning("Converting D to a dense matrix, because it is row rank deficient.")
      D = as.matrix(D)
      sparse = FALSE
    }
  }
  
  # If we got here, D is either tall or wide and row 
  # rank deficient
  
  if (sparse & (n < m)) {
    warning("Converting D to a dense matrix, because it has more rows than columns.")
    D = as.matrix(D)
    sparse = FALSE
  }

  # Compute a QR factorization of D
  x = qr(D)
  D = D[,x$pivot] # Pivot the columns of D
  y = y[x$pivot]  # Pivot the entries of y
  R = qr.R(x,complete=TRUE) # m x n
  
  i = which(abs(diag(R))<rtol)
  if (length(i)==0) k = 0
  else k = min(m,n)-min(i)+1
  
  # Do Givens rotations on the columns of R to make it
  # upper triangular. We also have to rotate y and D
  a = maketri2(y,D,R,k)
  y = a$y                   
  D = a$D                  
  R = a$R
  
  # This is the number of columns of R (and hence D),
  # from the left, that are zero
  q = k+max(n-m,0)

  # Trim y and D
  y = y[Seq(q+1,n)]                       # (n-q) x 1
  D = D[,Seq(q+1,n),drop=FALSE]           # m x (n-q)

  # Note: here it is necessarily true that m>n-q, i.e.,
  # that the resulting D is necessarily tall. We must have
  # m>=n-q because it has full column-rank, and we cannot
  # have m=n-q because it was row-rank deficient (otherwise
  # it would have been handled by the wide code, above)
  
  # Carry on with our current QR
  R = R[Seq(1,n-q),Seq(q+1,n),drop=FALSE] # (n-q) x (n-q)
  Q = qr.Q(x,complete=TRUE)               # We need the full Q
  Q1 = Q[,Seq(1,n-q),drop=FALSE]          # m x (n-q) 
  Q2 = Q[,Seq(n-q+1,m),drop=FALSE]        # m x (m-n+q)
    
  out = dualpathTall(y,D,Q1,Q2,R,q,approx,maxsteps,minlam,rtol,btol,verbose)
 
  # Construct beta, fit, y, bls, while accounting for the fact that
  # we may have had zero columns in D 
  out$df = out$df + coldif
  beta = matrix(y0,n0,length(out$lambda))
  beta[j,] = y0[j] - t(D0[,j])%*%out$u  
  colnames(beta) = colnames(out$u)
  out$beta = beta
  out$fit = beta
  out$y = y0
  out$bls = y0

  # Add to pathobjs component
  out$pathobjs$n0 = n0
  out$pathobjs$y0 = y0
  out$pathobjs$j = j
  out$pathobjs$D0 = D0
  out$pathobjs$coldif = coldif
  return(out)
}
