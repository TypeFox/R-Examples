dualpathFusedX.mod <- function(y, X, D, approx=FALSE, maxsteps=2000, minlam=0,
                            tol=1e-11) {
  m = nrow(D)
  p = ncol(D)
  n = length(y)
  
  
  # Find the minimum 2-norm solution, using some linear algebra 
  # tricks and a little bit of graph theory
  #save(D, file = 'D.RData')
  #temp = is.numeric(D[1,])
  #save(temp, file = 'temp.RData')
  #D = as.matrix(D)
  #L = abs(crossprod(D))
  L = abs(t(D) %*% D)
  diag(L) = 0
  gr = graph.adjacency(L,mode="upper") # Underlying graph
  cl = clusters(gr)                         
  q = cl$no                            # Number of clusters
  i = cl$membership                    # Cluster membership

  # First we project y onto the row space of D*X^+
  xy = t(X)%*%y
  A = matrix(0,n,q)
  z = numeric(q)
  
  # For efficiency, don't loop over singletons
  tab = tabulate(i)
  oo = which(tab[i]==1)
  if (length(oo)>0) {
    j = i[oo]
    A[,j] = X[,oo,drop=FALSE] 
    z[j] = xy[oo]
  }
  
  # Now consider all groups with at least two elements
  grps = which(tab>1)
  for (j in grps) {
    oo = which(i==j)
    A[,j] = rowMeans(X[,oo,drop=FALSE])
    z[j] = mean(xy[oo])
  }
  
  #e = ginv(crossprod(A)) %*% z
  e = ginv(t(A)%*%A) %*% z
  g = xy-t(X)%*%(A%*%e)

  # Here we perform our usual fused lasso solve but
  # with g in place of y
  x = f = numeric(p)
  
  # Again for efficiency, don't loop over singletons
  oo = which(tab[i]==1)
  if (length(oo)>0) {
    f[oo] = e[i][oo]
  }
  
  # Same for groups with two elements (doubletons?)
  oi = order(i)
  oo = which(tab[i][oi]==2)
  if (length(oo)>0) {
    f[oi][oo] = e[i][oi][oo]/2
    mm = colMeans(matrix(g[oi][oo],nrow=2))
    ii = oo[Seq(1,length(oo),by=2)]
    x[oi][ii] = g[oi][ii] - mm
  }
  
  # Now all groups with at least three elements
  cs = cumsum(tab)
  grps = which(tab>2)
  for (j in grps) {
    oo = oi[Seq(cs[j]-tab[j]+1,cs[j])]    
    f[oo] = e[j]/length(oo)
    gj = g[oo]
    Lj = t(D[,oo[-1]])%*% D[,oo[-1]]
    #Lj = crossprod(Matrix(D[,oo[-1]],sparse=TRUE))
    x[oo][-1] = as.numeric(ginv(as.matrix(Lj)) %*% (gj-mean(gj))[-1])
    
  }

  uhat = as.numeric(D%*%x)     # Dual solution
  betahat = f                  # Primal solution
  ihit = which.max(abs(uhat))  # Hitting coordinate
  hit = abs(uhat[ihit])        # Critical lambda
  s = sign(uhat[ihit])         # Sign

  # Now iteratively find the new dual solution, and
  # the next critical point
  
  # Things to keep track of, and return at the end
  buf = min(maxsteps,1000)
  lams = numeric(buf)        # Critical lambdas
  h = logical(buf)           # Hit or leave?
  df = numeric(buf)          # Degrees of freedom

  lams[1] = hit
  h[1] = TRUE
  df[1] = q
  
  # We only record the solutions if there is no
  # filebacking
  u = matrix(0,m,buf)      # Dual solutions
  beta = matrix(0,p,buf)   # Primal solutions
  u[,1] = uhat
  beta[,1] = betahat
  
  # Update our graph
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
  
  # Other things to keep track of, but not return
  r = 1                      # Size of boundary set
  B = ihit                   # Boundary set
  I = Seq(1,m)[-ihit]        # Interior set
  D1 = D[-ihit,,drop=FALSE]  # Matrix D[I,]
  D2 = D[ihit,,drop=FALSE]   # Matrix D[B,]
  k = 2                      # What step are we at?

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

    # For efficiency, don't loop over singletons
    tab = tabulate(i)
    oo = which(tab[i]==1)
    if (length(oo)>0) {
      j = i[oo]
      A[,j] = X[,oo,drop=FALSE] 
      z[j,1] = xy[oo]
      z[j,2] = Ds[oo]
    }
    
    # Now consider all groups with at least two elements
    grps = which(tab>1)
    for (j in grps) {
      oo = which(i==j)
      A[,j] = rowMeans(X[,oo,drop=FALSE])
      z[j,1] = mean(xy[oo])
      z[j,2] = mean(Ds[oo])
    }
    #e = ginv(crossprod(A)) %*% z
    e = ginv(t(A)%*%A) %*% z
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
      if (length(oo)>0) {
        fa[oo] = ea[i][oo]
        fb[oo] = eb[i][oo]
      }
      
      # Same for groups with two elements (doubletons?)
      oi = order(i)
      oo = which(tab[i][oi]==2)
      if (length(oo)>0) {
        fa[oi][oo] = ea[i][oi][oo]/2
        fb[oi][oo] = eb[i][oi][oo]/2
        ma = colMeans(matrix(ga[oi][oo],nrow=2))
        mb = colMeans(matrix(gb[oi][oo],nrow=2))
        ii = oo[Seq(1,length(oo),by=2)]
        xa[oi][ii] = ga[oi][ii] - ma
        xb[oi][ii] = gb[oi][ii] - mb
      }
      
      # Now all groups with at least three elements
      cs = cumsum(tab)
      grps = which(tab>2)
      for (j in grps) {
        oo = oi[Seq(cs[j]-tab[j]+1,cs[j])]    
        fa[oo] = ea[j]/length(oo)
        fb[oo] = eb[j]/length(oo)
        gaj = ga[oo]
        gbj = gb[oo]
        #Lj = crossprod(Matrix(D1[,oo[-1]],sparse=TRUE))
        Lj = t(D1[,oo[-1]]) %*% D1[,oo[-1]]
        xa[oo][-1] = as.numeric(ginv(as.matrix(Lj)) %*% (gaj-mean(gaj))[-1])
        xb[oo][-1] = as.numeric(ginv(as.matrix(Lj)) %*% (gbj-mean(gaj))[-1])
        
        
      }
      
      a = as.numeric(D1%*%xa)
      b = as.numeric(D1%*%xb)
      shits = Sign(a)
      hits = a/(b+shits);

      # Make sure none of the hitting times are larger
      # than the current lambda (precision issue)
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
      df[k] = q
      uhat = numeric(m)
      uhat[B] = hit*s
      uhat[I] = a-hit*b
      betahat = fa-hit*fb      

      # Update our graph
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
      
      # Update all other variables
      r = r+1
      B = c(B,I[ihit])
      I = I[-ihit]
      s = c(s,shit)
      D2 = rBind(D2,D1[ihit,])
      D1 = D1[-ihit,,drop=FALSE]
        
    }
              
    # Otherwise a leaving time comes next
    else {
      # Record the critical lambda and properties
      lams[k] = leave
      h[k] = FALSE
      df[k] = q
      uhat = numeric(m)
      uhat[B] = leave*s
      uhat[I] = a-leave*b
      betahat = fa-leave*fb
      
      # Update our graph
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
      
      # Update all other variables
      r = r-1
      I = c(I,B[ileave])
      B = B[-ileave]
      s = s[-ileave]
      D1 = rBind(D1,D2[ileave,])
      D2 = D2[-ileave,,drop=FALSE]

      
    }

    # Only record the solutions if we are not
    # filebacking
      u[,k] = uhat
      beta[,k] = betahat
    
    # Step counter
    k = k+1
  }

  # Trim 
  lams = lams[Seq(1,k-1)]
  h = h[Seq(1,k-1)]
  df = df[Seq(1,k-1)]
  u = u[,Seq(1,k-1),drop=FALSE]
  beta = beta[,Seq(1,k-1),drop=FALSE]

  # If we reached the maximum number of steps
  if (k>maxsteps) completepath = FALSE

  # If we reached the minimum lambda
  else if (lams[k-1]<minlam) completepath = FALSE

  # Otherwise, note that we completed the path
  else completepath = TRUE

  # The least squares solution (lambda=0)
  bls = NULL
  if (completepath) bls = fa
  
    colnames(u) = as.character(round(lams,3))
    colnames(beta) = as.character(round(lams,3))
    return(list(lambda=lams,beta=beta,fit=X%*%beta,u=u,hit=h,df=df,y=y,X=X,
                completepath=completepath,bls=bls,alpha=0))
}
