# We compute a solution path of the generalized lasso dual problem:
#
# \hat{u}(\lambda) =
# \argmin_u \|y - D^T u\|_2^2 \rm{s.t.} \|\u\|_\infty \leq \lambda
#
# where D is m x n. Here we treat the "tall" case, where m > n and
# rank(D) = n.
#
# Note: the df estimates at each lambda_k can be thought of as the df
# for all solutions corresponding to lambda in (lambda_k,lambda_{k-1}),
# the open interval to the *right* of the current lambda_k.

dualpathTall <- function(y, D, Q1, Q2, R, q0, approx=FALSE, maxsteps=2000,
                         minlam=0, rtol=1e-7, btol=1e-7, verbose=FALSE,
                         object=NULL) {
  # If we are starting a new path
  if (is.null(object)) {
    m = nrow(D)
    n = ncol(D)

    # We are given a QR factorization of D. The
    # dimensions are:
    # y:  n x 1
    # D:  m x n
    # Q1: m x n
    # Q2: m x (m-n)
    # R:  n x n

    # This is the number of columns of R (and hence D),
    # from the left, that are zero
    q = 0
    
    # Find the minimum 2-norm solution, and find the 
    # first critical point
    z = Backsolve(R,y,q)
    uhat = Q1%*%z                # Dual solution
    ihit = which.max(abs(uhat))  # Hitting coordinate
    hit = abs(uhat[ihit])        # Critical lambda
    s = sign(uhat[ihit])         # Sign

    if (verbose) {
      cat(sprintf("1. lambda=%.3f, adding coordinate %i, |B|=%i...",
                  hit,ihit,1))
    }

    # Now iteratively find the new dual solution, and
    # the next critical point

    # Things to keep track of, and return at the end
    buf = min(maxsteps,1000)
    u = matrix(0,m,buf)        # Dual solutions
    lams = numeric(buf)        # Critical lambdas
    h = logical(buf)           # Hit or leave?
    df = numeric(buf)          # Degrees of freedom

    lams[1] = hit
    h[1] = TRUE
    df[1] = q0
    u[,1] = uhat
    
    # Other things to keep track of, but not return
    r = 1                      # Size of boundary set
    B = ihit                   # Boundary set
    I = Seq(1,m)[-ihit]        # Interior set
    D1 = D[-ihit,,drop=FALSE]  # Matrix D[I,]
    D2 = D[ihit,,drop=FALSE]   # Matrix D[B,]
    k = 2                      # What step are we at?
    
    # Throughout the algorithm, D1 = Q1*R, and the
    # dimensions are:
    # D1: (m-r) x n
    # D2: r x n
    # Q1: (m-r) x min(n,m-r)
    # Q2: (m-r) x min(m-r-n,0)
    # R:  min(n,m-r) x n
    # Remember that the first q columns of D1 and R
    # are identically zero
  }

  # If iterating an already started path
  else { 
    # Grab variables needed to construct the path
    lambda = NULL
    for (j in 1:length(object)) {
      if (names(object)[j] != "pathobjs") {
        assign(names(object)[j], object[[j]])
      }
    }
    for (j in 1:length(object$pathobjs)) {
      assign(names(object$pathobjs)[j], object$pathobjs[[j]])
    }
    lams = lambda
  }

  tryCatch({  
    while (k<=maxsteps && lams[k-1]>=minlam) {
      ##########
      # Check if we've reached the end of the buffer
      if (k > length(lams)) {
        buf = length(lams)
        lams = c(lams,numeric(buf))
        h = c(h,logical(buf))
        df = c(df,numeric(buf))
        u = cbind(u,matrix(0,m,buf))
      }

      ##########
      # If we just added an element to B, downdate the
      # QR factorization
      if (h[k-1]) {
        x = downdateT(Q1,Q2,R,ihit)
        Q1 = x$Q1
        Q2 = x$Q2
        R = x$R
        
        # If D1 changed rank after we removed a row from
        # it, then we have to make R tridiagonal again. 
        # We need to consider three cases. When we talk
        # below about the dimension or diagonal of R, we
        # mean after having removed the zero padding

        # If R has no rows, then it changed rank
        if (nrow(R)==0) {
          # Increment the zero-column counter
          q = q+1
        }

        # If R has one less row than column, then it 
        # changed rank
        else if (nrow(R)==n-q-1) {
          # Do Givens rotations on the columns of R to make
          # it tridiagonal. We also have to rotate y,D1,D2
          x = maketri3(y,D1,D2,R,q)
          y = x$y                   
          D1 = x$D1
          D2 = x$D2
          R = x$R

          # Increment the zero-column counter
          q = q+1
        }

        # If R is square
        else if (nrow(R)>=n-q) {
          # If R has a zero on the diagonal, then it
          # changed rank
          d = diag(R[Seq(1,n-q),Seq(q+1,n),drop=FALSE])
          if (min(abs(d))<rtol) {
            # Do Givens rotations on the rows and columns
            # of R to make it tridiagonal. We also have to
            # rotate y,D1,D2 and Q1,Q2
            i = max(which(abs(d)<rtol))
            x = maketri4(y,D1,D2,Q1,Q2,R,q,i)
            y = x$y                   
            D1 = x$D1
            D2 = x$D2
            Q1 = x$Q1
            Q2 = x$Q2
            R = x$R
          
            # Increment the zero-column counter
            q = q+1
          }
        }
      }
      
      # Otherwise we just deleted an element from B, so
      # update the QR factorization
      else {
        x = updateT(y,D1,D2,Q1,Q2,R,q,D1[m-r,])
        y = x$y
        D1 = x$D1
        D2 = x$D2
        Q1 = x$Q1
        Q2 = x$Q2
        R = x$R

        # If D1 didn't change rank after we added a row 
        # to it, then we have to make R tridiagonal again.
        
        # If q is zero or the first element on the diagonal
        # of R (after having removed its first q columns) is
        # zero, then it didn't change rank
        if (q==0 || abs(R[1,q])<rtol) {
          # Do Givens rotations on the rows of R to make it
          # tridiagonal. We also have to rotate Q1,Q2 (and
          # pass y,D1,D2 even though they will not change)
          x = maketri4(y,D1,D2,Q1,Q2,R,q,1) 
          y = x$y                   
          D1 = x$D1
          D2 = x$D2
          Q1 = x$Q1
          Q2 = x$Q2
          R = x$R
        }

        # Otherwise, R changed rank, so decrement the zero-
        # column counter
        else {
          q = q-1
        }
      }
      
      # If the R factor is degenerate (it has a zero on 
      # the diagonal), then we just re-factor completely
      if (r<m && min(abs(diag(R[Seq(1,n-q),Seq(q+1,n),drop=FALSE])))<rtol) {
        if (verbose) {
          cat("\n(Recomputing QR factorization, for numerical stability...)")
        }
        
        x = qr(D1[,Seq(q+1,n)])
        Q = qr.Q(x,complete=TRUE)                    
        Q1 = Q[,Seq(1,min(n,m-r)),drop=FALSE]                          
        Q2 = Q[,Seq(min(n,m-r)+1,m-r),drop=FALSE]
        R = qr.R(x,complete=TRUE)
        R = cbind(matrix(0,min(n,m-r),q),R[Seq(1,min(n,m-r)),])
      }

      ##########
      # (Unlike in dualpathWide, we aren't storing this)
      Ds = t(D2)%*%s
      
      # If the interior is empty, then nothing will hit
      if (r==m) {
        a = b = numeric(0)
        hit = 0
      }
      
      # Otherwise, find the next hitting time
      else {
        z = Backsolve(R,y,q)
        a = Q1%*%z
        z = Backsolve(R,Ds,q)
        b = Q1%*%z
        shits = Sign(a)
        hits = a/(b+shits);

        # Make sure none of the hitting times are larger
        # than the current lambda (precision issue)
        hits[hits>lams[k-1]+btol] = 0
        hits[hits>lams[k-1]] = lams[k-1]
        
        ihit = which.max(hits)
        hit = hits[ihit]
        shit = shits[ihit]
      }
      
      ##########
      # If nothing is on the boundary, then nothing will leave
      # Also, skip this if we are in "approx" mode
      if (r==0 || approx) {
        leave = 0
      }
      
      # Otherwise, find the next leaving time
      else {
        # Note: c and d involve projecting a vector onto the null
        # space of D_{-B}. We can run into trouble when a vector
        # is orthogonal to this null space, so these should be zero,
        # but due to numerical inaccuracy it's simply very close to
        # zero. dualpathWide doesn't suffer from this problem, because
        # in its case D_{-B} is always full row rank so the null space
        # is trivial. And dualpathFused doesn't suffer from this problem
        # because we can do the projection explicitly, i.e. by averaging 
        # over the connected components of the underlying graph. 
        c = s*(D2%*%(y-t(D1)%*%a))
        d = s*(D2%*%(Ds-t(D1)%*%b))
        leaves = c/d
        
        # c must be negative 
        leaves[c>=0] = 0
        
        # Make sure none of the leaving times are larger
        # than the current lambda (precision issue)
        leaves[leaves>lams[k-1]+btol] = 0 
        leaves[leaves>lams[k-1]] = lams[k-1]
        
        ileave = which.max(leaves)
        leave = leaves[ileave]
      }

      ##########
      # Stop if the next critical point is negative
      if (hit<=0 && leave<=0) break

      # If a hitting time comes next
      if (hit > leave) {
        # Record the critical lambda and solution
        lams[k] = hit
        h[k] = TRUE
        df[k] = q0+q
        uhat = numeric(m)
        uhat[B] = hit*s
        uhat[I] = a-hit*b
        u[,k] = uhat
        
        # Update all of the variables
        r = r+1
        B = c(B,I[ihit])
        I = I[-ihit]
        s = c(s,shit)
        D2 = rbind(D2,D1[ihit,])
        D1 = D1[-ihit,,drop=FALSE]
          
        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, adding coordinate %i, |B|=%i...",
                      k,hit,B[r],r))
        }
      }
                
      # Otherwise a leaving time comes next
      else {
        # Record the critical lambda and solution
        lams[k] = leave
        h[k] = FALSE
        df[k] = q0+q
        uhat = numeric(m)
        uhat[B] = leave*s
        uhat[I] = a-leave*b
        u[,k] = uhat
        
        # Update all of the variables
        r = r-1
        I = c(I,B[ileave])
        B = B[-ileave]
        s = s[-ileave]
        D1 = rbind(D1,D2[ileave,])
        D2 = D2[-ileave,,drop=FALSE]

        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, deleting coordinate %i, |B|=%i...",
                      k,leave,I[m-r],r))
        }
      }

      # Step counter
      k = k+1
    }
  }, error = function(err) {
    err$message = paste(err$message,"\n(Path computation has been terminated;",
      " partial path is being returned.)",sep="")
    warning(err)})

  # Trim
  lams = lams[Seq(1,k-1)]
  h = h[Seq(1,k-1)]
  df = df[Seq(1,k-1),drop=FALSE]
  u = u[,Seq(1,k-1),drop=FALSE]

  # Save needed elements for continuing the path
  pathobjs = list(type="tall", r=r, B=B, I=I, Q1=Q1, approx=approx,
    Q2=Q2, k=k, df=df, D1=D1, D2=D2, ihit=ihit, m=m, n=n, q=q, h=h,
    R=R, q0=q0, rtol=rtol, btol=btol, s=s, y=y)
  
  # If we reached the maximum number of steps
  if (k>maxsteps) {
    if (verbose) {
      cat(sprintf("\nReached the maximum number of steps (%i),",maxsteps))
      cat(" skipping the rest of the path.")
    }
    completepath = FALSE
  }

  # If we reached the minimum lambda
  else if (lams[k-1]<minlam) {
    if (verbose) {
      cat(sprintf("\nReached the minimum lambda (%.3f),",minlam))
      cat(" skipping the rest of the path.")
    }
    completepath = FALSE
  }
  
  # Otherwise, note that we completed the path
  else completepath = TRUE
  
  if (verbose) cat("\n")

  # The parent funtion will return the proper beta, fit, y, bls
  colnames(u) = as.character(round(lams,3))
  
  return(list(lambda=lams,beta=NA,fit=NA,u=u,hit=h,df=df,y=NA,
              completepath=completepath,bls=NA,pathobjs=pathobjs))
}
