# We compute a solution path of the sparse fused lasso dual problem:
#
# \hat{u}(\lambda) =
# \argmin_u \|y - (X^+)^T D^T u\|_2^2 \rm{s.t.} \|\u\|_\infty \leq \lambda
#
# where D is (a multiple of) incidence matrix of a given graph, row-
# binded with (a multiple of) the identity matrix, and X is a full column
# rank predictor matrix, X^+ being its pseudoinverse.
#
# Fortuitously, we never have to fully invert X (i.e. compute its pseudo-
# inverse).
#
# Note: the df estimates at each lambda_k can be thought of as the df
# for all solutions corresponding to lambda in (lambda_k,lambda_{k-1}),
# the open interval to the *right* of the current lambda_k.

dualpathFusedL1X <- function(y, X, D, D0, gamma, approx=FALSE, maxsteps=2000,
                             minlam=0, rtol=1e-7, btol=1e-7, eps=1e-4,
                             verbose=FALSE, object=NULL) {
  # If we are starting a new path
  if (is.null(object)) {
    m = nrow(D)
    p = ncol(D)
    n = length(y)
    numedges = m-p
    numnodes = p

    # Modify y,X,n in the case of a ridge penalty, but
    # keep the originals
    y0 = y
    X0 = X
    if (eps>0) {
      y = c(y,rep(0,p))
      X = rbind(X,diag(sqrt(eps),p))
      n = n+p
    }

    # Find the minimum 2-norm solution, using some linear algebra
    # tricks and a little bit of graph theory
    L = abs(crossprod(D0))
    diag(L) = 0
    gr = graph.adjacency(L,mode="upper") # Underlying graph
    cl = clusters(gr)
    q = cl$no                            # Number of clusters
    i = cl$membership                    # Cluster membership

    # First we project y onto the row space of D*X^+
    xy = t(X)%*%y
    g = xy

    # Here we perform our usual fused lasso solve but
    # with g in place of y
    x = numeric(p)

    # For efficiency, don't loop over singletons
    tab = tabulate(i)
    oo = which(tab[i]==1)
    if (length(oo)>0) {
      x[oo] = g[oo]/(diag(L)[oo])
      }

    # Now all groups with at least two elements
    oi = order(i)
    cs = cumsum(tab)
    grps = which(tab>1)
    for (j in grps) {
      oo = oi[Seq(cs[j]-tab[j]+1,cs[j])]
      Lj = crossprod(Matrix(D[,oo],sparse=TRUE))
      x[oo] = as.numeric(solve(Lj,g[oo]))
    }

    uhat = as.numeric(D%*%x)     # Dual solution
    betahat = numeric(p)         # Primal solution
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

    lams[1] = hit
    h[1] = TRUE
    df[1] = 0

    u = matrix(0,m,buf)      # Dual solutions
    beta = matrix(0,p,buf)   # Primal solutions
    u[,1] = uhat
    beta[,1] = betahat

    # Special interior set over nodes
    I0 = rep(TRUE,numnodes)

    # Update the graph if we need to, otherwise
    # update the special interior set
    if (ihit <= numedges) {
      ed = which(D[ihit,]!=0)
      gr[ed[1],ed[2]] = 0             # Delete edge
      newcl = subcomponent(gr,ed[1])  # New cluster
      oldcl = which(i==i[ed[1]])      # Old cluster
      # If these two clusters aren't the same, update
      # the memberships
      if (length(newcl)!=length(oldcl) || any(sort(newcl)!=sort(oldcl))) {
        i[newcl] = q+1
        q = q+1
      }
    }
    else {
      I0[ihit-numedges] = FALSE
    }

    # Other things to keep track of, but not return
    r = 1                      # Size of boundary set
    B = ihit                   # Boundary set
    I = Seq(1,m)[-ihit]        # Interior set
    D1 = D[-ihit,,drop=FALSE]  # Matrix D[I,]
    D2 = D[ihit,,drop=FALSE]   # Matrix D[B,]
    k = 2                      # What step are we at?
  }

  # If iterating already started path
  else {
    # Grab variables from outer object
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
      Ds = as.numeric(t(D2)%*%s)

      # Precomputation for the hitting times: first we project
      # y and Ds onto the row space of D1*X^+
      A = matrix(0,n,q)
      z = matrix(0,q,2)
      nz = rep(FALSE,q)

      # For efficiency, don't loop over singletons
      tab = tabulate(i)
      oo = which(tab[i]==1)
      oo2 = oo[!I0[oo]]
      if (length(oo2)>0) {
        j = i[oo2]
        A[,j] = X[,oo2,drop=FALSE]
        z[j,1] = xy[oo2]
        z[j,2] = Ds[oo2]
        nz[j] = TRUE
      }

      # Now consider all groups with at least two elements
      grps = which(tab>1)
      for (j in grps) {
        oo = which(i==j)
        if (all(!I0[oo])) {
          A[,j] = rowMeans(X[,oo,drop=FALSE])
          z[j,1] = mean(xy[oo])
          z[j,2] = mean(Ds[oo])
          nz[j] = TRUE
        }
      }

      nzq = sum(nz)
      e = matrix(0,q,2)
      if (nzq>0) {
        R = qr.R(qr(A[,nz]))
        e[nz,] = backsolve(R,forwardsolve(R,z[nz,,drop=FALSE],upper.tri=TRUE,transpose=TRUE))
        # Note: using a QR here is preferable than simply calling
        # e[nz,] = solve(crossprod(A[,nz]),z[nz,,drop=FALSE]), for
        # numerical stablity. Plus, it's not really any slower
      }
      ea = e[,1]
      eb = e[,2]
      ga = xy-t(X)%*%(A%*%ea)
      gb = Ds-t(X)%*%(A%*%eb)

      # If the interior is empty, then nothing will hit
      if (r==m) {
        fa = ea[i]
        fb = eb[i]
        hit = 0
      }

      # Otherwise, find the next hitting time
      else {
        # Here we perform our usual fused lasso solve but
        # with ga in place of y and gb in place of Ds
        xa = xb = numeric(p)
        fa = fb = numeric(p)

        # For efficiency, don't loop over singletons
        oo = which(tab[i]==1)
        fa[oo] = ea[i][oo]
        fb[oo] = eb[i][oo]
        oo1 = oo[I0[oo]]
        if (length(oo1)>0) {
          Ldiag = diag(crossprod(Matrix(D1[,oo1],sparse=TRUE)))
          xa[oo1] = ga[oo1]/Ldiag
          xb[oo1] = gb[oo1]/Ldiag
        }

        # Now all groups with at least two elements
        oi = order(i)
        cs = cumsum(tab)
        grps = which(tab>1)
        for (j in grps) {
          oo = oi[Seq(cs[j]-tab[j]+1,cs[j])]
          fa[oo] = ea[j]/length(oo)
          fb[oo] = eb[j]/length(oo)
          gaj = ga[oo]
          gbj = gb[oo]

          if (any(I0[oo])) {
            Lj = crossprod(Matrix(D1[,oo],sparse=TRUE))
            xa[oo] = as.numeric(solve(Lj,gaj))
            xb[oo] = as.numeric(solve(Lj,gbj))
          }
          else {
            Lj = crossprod(Matrix(D1[,oo[-1]],sparse=TRUE))
            xa[oo][-1] = as.numeric(solve(Lj,(gaj-mean(gaj))[-1]))
            xb[oo][-1] = as.numeric(solve(Lj,(gbj-mean(gbj))[-1]))
          }
        }

        a = as.numeric(D1%*%xa)
        b = as.numeric(D1%*%xb)
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
        df[k] = nzq
        uhat = numeric(m)
        uhat[B] = hit*s
        uhat[I] = a-hit*b
        betahat = fa-hit*fb

        # Update our graph if we need to, otherwise
        # update the special interior set
        if (I[ihit] <= numedges) {
          ed = which(D1[ihit,]!=0)
          gr[ed[1],ed[2]] = 0             # Delete edge
          newcl = subcomponent(gr,ed[1])  # New cluster
          oldcl = which(i==i[ed[1]])      # Old cluster
          # If these two clusters aren't the same, update
          # the memberships
          if (length(newcl)!=length(oldcl) || any(sort(newcl)!=sort(oldcl))) {
            i[newcl] = q+1
            q = q+1
          }
        }
        else {
          I0[I[ihit]-numedges] = FALSE
        }

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
        df[k] = nzq
        uhat = numeric(m)
        uhat[B] = leave*s
        uhat[I] = a-leave*b
        betahat = fa-leave*fb

        # Update our graph if we need to, otherwise
        # update the special interior set
        if (B[ileave] <= numedges) {
          ed = which(D2[ileave,]!=0)
          gr[ed[1],ed[2]] = 1             # Add edge
          newcl = subcomponent(gr,ed[1])  # New cluster
          oldcl = which(i==i[ed[1]])      # Old cluster
          # If these two clusters aren't the same, update
          # the memberships
          if (length(newcl)!=length(oldcl) || !all(sort(newcl)==sort(oldcl))) {
            newno = i[ed[2]]
            oldno = i[ed[1]]
            i[oldcl] = newno
            i[i>oldno] = i[i>oldno]-1
            q = q-1
          }
        }
        else {
          I0[B[ileave]-numedges] = TRUE
        }

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

  # Save needed elements for continuing the path
  pathobjs = list(type="fused.l1.x" ,r=r, B=B, I=I, Q1=NA, approx=approx,
    Q2=NA, k=k, df=df, D1=D1, D2=D2, Ds=Ds, ihit=ihit, m=m, n=n, p=p, q=q,
    h=h, q0=NA, rtol=rtol, btol=btol, eps=eps, s=s, y=y, gr=gr, i=i,
    numedges=numedges, I0=I0, xy=xy)

  colnames(u) = as.character(round(lams,3))
  colnames(beta) = as.character(round(lams,3))
  return(list(lambda=lams,beta=beta,fit=X0%*%beta,u=u,hit=h,df=df,y=y0,X=X0,
              completepath=completepath,bls=bls,gamma=gamma,pathobjs=pathobjs))
}
