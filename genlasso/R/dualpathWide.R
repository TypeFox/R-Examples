# We compute a solution path of the generalized lasso dual problem:
#
# \hat{u}(\lambda) =
# \argmin_u \|y - D^T u\|_2^2 \rm{s.t.} \|\u\|_\infty \leq \lambda
#
# where D is m x n. Here we treat the "wide" case, where m <= n and
# rank(D) = m.
#
# Note: the df estimates at each lambda_k can be thought of as the df
# for all solutions corresponding to lambda in (lambda_k,lambda_{k-1}),
# the open interval to the *right* of the current lambda_k.

dualpathWide <- function(y, D, Q1, Q2, R, approx=FALSE, maxsteps=2000,
                         minlam=0, rtol=1e-7, btol=1e-7, verbose=FALSE,
                         object=NULL) {
  # If we are starting a new path
  if (is.null(object)) {
    m = nrow(D)
    n = ncol(D)

    # We are given a QR factorization of D^T. The
    # dimensions are:
    # y:  n x 1
    # D:  m x n
    # Q1: n x m
    # Q2: n x 0
    # R:  m x m
    
    # Compute the dual solution at infinity, and
    # find the first critical point
    uhat = backsolve(R,t(Q1)%*%y) # Dual solution
    ihit = which.max(abs(uhat))   # Hitting coordinate
    hit = abs(uhat[ihit])         # Critical lambda
    s = Sign(uhat[ihit])          # Sign
    
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
    df[1] = n-m
    u[,1] = uhat
    
    # Other things to keep track of, but not return
    r = 1                      # Size of boundary set
    B = ihit                   # Boundary set
    I = Seq(1,m)[-ihit]        # Interior set
    Ds = D[ihit,]*s            # Vector t(D[B,])%*%s
    D1 = D[-ihit,,drop=FALSE]  # Matrix D[I,]
    D2 = D[ihit,,drop=FALSE]   # Matrix D[B,]
    k = 2                      # What step are we at?
    
    # Throughout the algorithm, D1^T = Q1*R, and the 
    # dimensions are:
    # D1: n x (m-r)
    # D2: n x r
    # Q1: n x (m-r)
    # Q2: n x r
    # R:  (m-r) x (m-r)
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
          x = downdateW(Q1,Q2,R,ihit)
      } else {
      # Otherwise we just deleted an element from B, so
      # update the QR factorization
       x = updateW(Q1,Q2,R,D1[m-r,])
      }

      Q1 = x$Q1
      Q2 = x$Q2
      R = x$R
      
      # If the R factor is degenerate (it has a zero on 
      # the diagonal), then we just re-factor completely
      if (r<m && min(abs(diag(R)))<rtol) {
        x = qr(t(rbind(D1,D2)))
        Q = qr.Q(x,complete=FALSE)                    
        Q1 = Q[,Seq(1,m-r),drop=FALSE]                          
        Q2 = Q[,Seq(m-r+1,m),drop=FALSE]              
        R = qr.R(x)[Seq(1,m-r),Seq(1,m-r),drop=FALSE] 
        
        if (verbose) {  
          cat("\n(Recomputing QR factorization, for numerical stability...)")
        }
      }
      
      ##########
      # If the interior is empty, then nothing will hit
      if (r==m) {
        a = b = numeric(0)
        hit = 0    
      }
    
      # Otherwise, find the next hitting time
      else {
        a = backsolve(R,t(Q1)%*%y) 
        b = backsolve(R,t(Q1)%*%Ds)
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
        df[k] = n-m+r
        uhat = numeric(m)
        uhat[B] = hit*s
        uhat[I] = a-hit*b
        u[,k] = uhat
        
        # Update all of the variables
        r = r+1
        B = c(B,I[ihit])
        I = I[-ihit]
        Ds = Ds + D1[ihit,]*shit
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
        df[k] = n-m+r
        uhat = numeric(m)
        uhat[B] = leave*s
        uhat[I] = a-leave*b
        u[,k] = uhat
        
        # Update all of the variables
        r = r-1
        I = c(I,B[ileave])
        B = B[-ileave]
        Ds = Ds - D2[ileave,]*s[ileave]
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
  pathobjs = list(type="wide", r=r, B=B, I=I, Q1=Q1, approx=approx,
    Q2=Q2, k=k, df=df, D1=D1, D2=D2, Ds=Ds, ihit=ihit, m=m, n=n, q=q,
    h=h, R=R, q0=NA, rtol=rtol, btol=btol, s=s, y=y)
  
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
