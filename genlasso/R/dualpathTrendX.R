# We compute a solution path of the trend filtering dual problem:
#
# \hat{u}(\lambda) =
# \argmin_u \|y - (X^+)^T D^T u\|_2^2 \rm{s.t.} \|\u\|_\infty \leq \lambda
#
# where D is a trend filtering matrix, and X has full column rank, X^+ being
# its pseudoinverse.
#
# Fortuitously, we never have to fully invert X (i.e. compute its pseudo-
# inverse).
#
# Note: the df estimates at each lambda_k can be thought of as the df
# for all solutions corresponding to lambda in (lambda_k,lambda_{k-1}),
# the open interval to the *right* of the current lambda_k.

dualpathTrendX <- function(y, pos, X, D, ord, approx=FALSE, maxsteps=2000,
                           minlam=0, rtol=1e-7, btol=1e-7, eps=1e-4,
                           verbose=FALSE, object=NULL) {
  # If we are starting a new path
  if (is.null(object)) {
    m = nrow(D)
    p = ncol(D)
    n = length(y)

    # Modify y,X in the case of a ridge penalty, but
    # keep the originals 
    y0 = y
    X0 = X
    if (eps>0) {
      y = c(y,rep(0,p))
      X = rbind(X,diag(sqrt(eps),p))
      n = n+p
    }
  
    # Find the minimum 2-norm solution, using some linear algebra 
    # tricks and a little bit of discrete calculus
    if (is.null(pos)) pos = 1:p
    Pos = matrix(rep(pos,each=ord),ord,p)
    basis = matrix(0,p,ord+1)    # Basis for null(D)
    basis[,1] = rep(1,p)
    for (i in Seq(2,ord+1)) {
      ii = Seq(1,i-1)
      basis[,i] = apply(pmax(Pos[ii,,drop=FALSE]-pos[ii],0),
             2,prod)/factorial(i-1) 
    }
    
    # First project onto the row space of D*X^+
    xy = t(X)%*%y
    A = X%*%basis
    z = t(basis)%*%xy
    R = qr.R(qr(A))
    e = backsolve(R,forwardsolve(R,z,upper.tri=TRUE,transpose=TRUE))
    # Note: using a QR here is preferable than simply calling
    # e = solve(crossprod(A),z), for numerical stablity. Plus,
    # it's not really any slower
    g = xy-t(X)%*%(A%*%e)
  
   # Here we perform our usual trend filter solve but
    # with g in place of y
    x = qr(t(D))
    uhat = backsolveSparse(x,g)  # Dual solution
    betahat = basis%*%e          # Primal solution
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
    lams = numeric(buf)        # Critical lambdas
    h = logical(buf)           # Hit or leave?
    df = numeric(buf)          # Degrees of freedom
    u = matrix(0,m,buf)        # Dual solutions
    beta = matrix(0,p,buf)     # Primal solutions
    
    lams[1] = hit
    h[1] = TRUE
    df[1] = ncol(basis)
    u[,1] = uhat
    beta[,1] = betahat

    # Update our basis
    newbv = apply(pmax(Pos-pos[Seq(ihit+1,ihit+ord)],0),
      2,prod)/factorial(ord)
    newbv[Seq(1,ihit+ord)] = 0 # Only needed when ord=0
    basis = cbind(basis,newbv)
  
    # Other things to keep track of, but not return
    r = 1                      # Size of boundary set
    B = ihit                   # Boundary set
    I = Seq(1,m)[-ihit]        # Interior set
    D1 = D[-ihit,,drop=FALSE]  # Matrix D[I,]
    D2 = D[ihit,,drop=FALSE]   # Matrix D[B,]
    k = 2                      # What step are we at?
  }

  # If iterating an already existing path
  else {
    # Grab variables from outer object to construct the path
    lambda = NULL
    for (j in 1:length(object)) {
      if (names(object)[j] != "pathobjs") {
        assign(names(object)[j], object[[j]])
      }
    }

    # Trick: save y,X from outer object
    y0 = y
    X0 = X

    # Grab variables from inner object
    for (j in 1:length(object$pathobjs)) {
      assign(names(object$pathobjs)[j], object$pathobjs[[j]])
    }
    lams = lambda

    # In the case of a ridge penalty, modify X
    if (eps>0) X = rbind(X,diag(sqrt(eps),p))
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
        beta = cbind(beta,matrix(0,p,buf))
      }

      ##########
      # No updating, just recompute these every time
      x = qr(t(D1))
      Ds = as.numeric(t(D2)%*%s)
      
      # Precomputation for the hitting times: first we project
      # y and Ds onto the row space of D1*X^+
      A = X%*%basis
      z = t(basis)%*%cbind(xy,Ds)
      R = qr.R(qr(A))
      e = backsolve(R,forwardsolve(R,z,upper.tri=TRUE,transpose=TRUE))
      # Note: using a QR here is preferable than simply calling
      # e = solve(crossprod(A),z), for numerical stablity. Plus,
      # it's not really any slower
      ea = e[,1]
      eb = e[,2]
      ga = xy-t(X)%*%(A%*%ea)
      gb = Ds-t(X)%*%(A%*%eb)
      fa = basis%*%ea
      fb = basis%*%eb
      
      # If the interior is empty, then nothing will hit
      if (r==m) {
        hit = 0
      }
      
      # Otherwise, find the next hitting time
      else {
        # Here we perform our usual trend filter solve but
        # with ga in place of y and gb in place of Ds
        a = backsolveSparse(x,ga)
        b = backsolveSparse(x,gb)
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
        c = as.numeric(s*(D2%*%fa))
        d = as.numeric(s*(D2%*%fb))
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
        # Record the critical lambda and properties
        lams[k] = hit
        h[k] = TRUE
        df[k] = ncol(basis)
        uhat = numeric(m)
        uhat[B] = hit*s
        uhat[I] = a-hit*b
        betahat = fa-hit*fb

        # Update our basis
        newbv = apply(pmax(Pos-pos[Seq(I[ihit]+1,I[ihit]+ord)],0),
          2,prod)/factorial(ord)
        newbv[Seq(1,I[ihit]+ord)] = 0 # Only needed when ord=0
        basis = cbind(basis,newbv)
        
        # Update all other variables
        r = r+1
        B = c(B,I[ihit])
        I = I[-ihit]
        s = c(s,shit)
        D2 = rBind(D2,D1[ihit,])
        D1 = D1[-ihit,,drop=FALSE]
          
        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, adding coordinate %i, |B|=%i...",
                      k,hit,B[r],r))
        }
      }
                
      # Otherwise a leaving time comes next
      else {
        # Record the critical lambda and properties
        lams[k] = leave
        h[k] = FALSE
        df[k] = ncol(basis)
        uhat = numeric(m)
        uhat[B] = leave*s
        uhat[I] = a-leave*b
        betahat = fa-leave*fb
        
        # Update our basis
        basis = basis[,-(ord+1+ileave)]
        
        # Update all other variables
        r = r-1
        I = c(I,B[ileave])
        B = B[-ileave]
        s = s[-ileave]
        D1 = rBind(D1,D2[ileave,])
        D2 = D2[-ileave,,drop=FALSE]

        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, deleting coordinate %i, |B|=%i...",
                      k,leave,I[m-r],r))
        }
      }

      u[,k] = uhat
      beta[,k] = betahat

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
  df = df[Seq(1,k-1)]
  u = u[,Seq(1,k-1),drop=FALSE]
  beta = beta[,Seq(1,k-1),drop=FALSE]
  
  # If we reached the maximum number of steps
  if (k>maxsteps) {
    if (verbose) {
      cat(sprintf("\nReached the max number of steps (%i),",maxsteps))
      cat(" skipping the rest of the path.")
    }
    completepath = FALSE
  }

  # If we reached the minimum lambda
  else if (lams[k-1]<minlam) {
    if (verbose) {
      cat(sprintf("\nReached the min lambda (%.3f),",minlam))
      cat(" skipping the rest of the path.")
    }
    completepath = FALSE
  }

  # Otherwise, note that we completed the path
  else completepath = TRUE

  # The least squares solution (lambda=0)
  bls = NULL
  if (completepath) bls = fa
  
  if (verbose) cat("\n")
  # Save elements needed for continuing the path
  pathobjs = list(type="trend.x", r=r, B=B, I=I, Q1=NA, approx=approx,
    Q2=NA, k=k, df=df, D1=D1, D2=D2, Ds=Ds, ihit=ihit, m=m, n=n, p=p,
    q=q, h=h, q0=NA, rtol=rtol, btol=btol, eps=eps, s=s, y=y, ord=ord,
    pos=pos, Pos=Pos, basis=basis, xy=xy)
  
  colnames(u) = as.character(round(lams,3))
  colnames(beta) = as.character(round(lams,3))
  return(list(lambda=lams,beta=beta,fit=X0%*%beta,u=u,hit=h,df=df,y=y0,X=X0,
              completepath=completepath,bls=bls,pathobjs=pathobjs))
}
