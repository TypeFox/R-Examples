# We compute a solution path of the fused lasso dual problem:
#
# \hat{u}(\lambda) =
# \argmin_u \|y - D^T u\|_2^2 \rm{s.t.} \|\u\|_\infty \leq \lambda
#
# where D is (a multiple of) incidence matrix of a given graph, row-
# binded with (a multiple of) the identity matrix.
#
# Note: the df estimates at each lambda_k can be thought of as the df
# for all solutions corresponding to lambda in (lambda_k,lambda_{k-1}),
# the open interval to the *right* of the current lambda_k.

dualpathFusedL1 <- function(y, D, D0, gamma, approx=FALSE, maxsteps=2000,
                            minlam=0, rtol=1e-7, btol=1e-7, verbose=FALSE,
                            object=NULL) {
  # If we are starting a new path
  if (is.null(object)) {
    m = nrow(D)
    n = ncol(D)
    numedges = m-n
    numnodes = n

    # Find the minimum 2-norm solution, using some linear algebra
    # tricks and a little bit of graph theory
    L = abs(crossprod(D0))
    diag(L) = 0
    gr = graph.adjacency(L,mode="upper") # Underlying graph
    cl = clusters(gr)
    q = cl$no                            # Number of clusters
    i = cl$membership                    # Cluster membership
    x = numeric(n)

    # For efficiency, don't loop over singletons
    tab = tabulate(i)
    oo = which(tab[i]==1)
    if (length(oo)>0) {
      x[oo] = y[oo]/(diag(L)[oo])
    }

    # Now all groups with at least two elements
    oi = order(i)
    cs = cumsum(tab)
    grps = which(tab>1)
    for (j in grps) {
      oo = oi[Seq(cs[j]-tab[j]+1,cs[j])]
      Lj = crossprod(Matrix(D[,oo],sparse=TRUE))
      x[oo] = as.numeric(solve(Lj,y[oo]))
    }

    uhat = as.numeric(D%*%x)     # Dual solution
    betahat = numeric(n)         # Primal solution
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
    beta = matrix(0,n,buf)   # Primal solutions
    u[,1] = uhat
    beta[,1] = betahat

    # Special interior set over nodes
    I0 = rep(TRUE,numnodes)

    # Update the graph if we need to, otherwise
    # update the special interior set
    if (ihit <= numedges) {
      e = which(D[ihit,]!=0)
      gr[e[1],e[2]] = 0             # Delete edge
      newcl = subcomponent(gr,e[1]) # New cluster
      oldcl = which(i==i[e[1]])     # Old cluster
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

    # Other things to keep track of
    r = 1                      # Size of boundary set
    B = ihit                   # Boundary set
    I = Seq(1,m)[-ihit]        # Interior set
    D1 = D[-ihit,,drop=FALSE]  # Matrix D[I,]
    D2 = D[ihit,,drop=FALSE]   # Matrix D[B,]
    k = 2                      # What step are we at?
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
        beta = cbind(beta,matrix(0,n,buf))
      }

      ##########
      Ds = as.numeric(t(D2)%*%s)
      nzq = q

      # If the interior is empty, then nothing will hit
      if (r==m) {
        fa = y
        fb = Ds
        hit = 0
      }

      # Otherwise, find the next hitting time
      else {
        xa = xb = numeric(n)
        fa = fb = numeric(n)

        # For efficiency, don't loop over singletons
        tab = tabulate(i)
        oo = which(tab[i]==1)
        oo1 = oo[I0[oo]]
        oo2 = oo[!I0[oo]]
        if (length(oo1)>0) {
          Ldiag = diag(crossprod(Matrix(D1[,oo1],sparse=TRUE)))
          xa[oo1] = y[oo1]/Ldiag
          xb[oo1] = Ds[oo1]/Ldiag
          nzq = nzq-length(oo1)
        }
        if (length(oo2)>0) {
          fa[oo2] = y[oo2]
          fb[oo2] = Ds[oo2]
        }

        # Now all groups with at least two elements
        oi = order(i)
        cs = cumsum(tab)
        grps = which(tab>1)
        for (j in grps) {
          oo = oi[Seq(cs[j]-tab[j]+1,cs[j])]
          yj = y[oo]
          Dsj = Ds[oo]

          if (any(I0[oo])) {
            Lj = crossprod(Matrix(D1[,oo],sparse=TRUE))
            xa[oo] = as.numeric(solve(Lj,yj))
            xb[oo] = as.numeric(solve(Lj,Dsj))
            nzq = nzq-1
          }
          else {
            fa[oo] = mean(yj)
            fb[oo] = mean(Dsj)
            Lj = crossprod(Matrix(D1[,oo[-1]],sparse=TRUE))
            xa[oo][-1] = as.numeric(solve(Lj,(yj-mean(yj))[-1]))
            xb[oo][-1] = as.numeric(solve(Lj,(Dsj-mean(Dsj))[-1]))
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
          e = which(D1[ihit,]!=0)
          gr[e[1],e[2]] = 0             # Delete edge
          newcl = subcomponent(gr,e[1]) # New cluster
          oldcl = which(i==i[e[1]])     # Old cluster
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
          e = which(D2[ileave,]!=0)
          gr[e[1],e[2]] = 1             # Add edge
          newcl = subcomponent(gr,e[1]) # New cluster
          oldcl = which(i==i[e[1]])     # Old cluster
          # If these two clusters aren't the same, update
          # the memberships
          if (length(newcl)!=length(oldcl) || !all(sort(newcl)==sort(oldcl))) {
            newno = i[e[2]]
            oldno = i[e[1]]
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

  # The least squares solution (lambda=0)
  bls = y
  if (verbose) cat("\n")

  # Save needed elements for continuing the path
  pathobjs = list(type="fused.l1",r=r, B=B, I=I, Q1=NA, approx=approx,
    Q2=NA, k=k, df=df, D1=D1, D2=D2, Ds=Ds, ihit=ihit, m=m, n=n, q=q,
    h=h, q0=NA, rtol=rtol, btol=btol, s=s, y=y, gr=gr, i=i,
    numedges=numedges, I0=I0)

  colnames(u) = as.character(round(lams,3))
  colnames(beta) = as.character(round(lams,3))
  return(list(lambda=lams,beta=beta,fit=beta,u=u,hit=h,df=df,y=y,
              completepath=completepath,bls=bls,gamma=gamma,pathobjs=pathobjs))
}
