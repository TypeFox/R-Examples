### main functions
# single change-point
gseg1 = function(n, E, n0=0.05*n, n1=0.95*n, pval.appr=TRUE, skew.corr=TRUE, pval.perm=FALSE, B=100){
  n0 = ceiling(n0)
  n1 = floor(n1)
  Ebynode = vector("list", n)  
  for(i in 1:n) Ebynode[[i]]=rep(0,0)
  for(i in 1:nrow(E)){
    Ebynode[[E[i,1]]] = c(Ebynode[[E[i,1]]],E[i,2])
    Ebynode[[E[i,2]]] = c(Ebynode[[E[i,2]]],E[i,1])
  }
  r1 = gcp1bynode(n,Ebynode,n0,n1)
  cat("Estimated change-point location:", r1$tauhat, "\n")
  cat("Zmax:", r1$Zmax, "\n")
  if (pval.appr==TRUE){
    mypval1 = pval1(n,E,Ebynode,r1$Zmax,skew.corr,n0,n1)
    r1$pval.appr = min(mypval1,1)
    cat("Approximated p-value:", r1$pval.appr, "\n")
  }
  if (pval.perm==TRUE){
    mypval2 = permpval1(n,Ebynode,r1$Zmax,B,n0,n1)
    r1$pval.perm = min(mypval2$pval,1)
    r1$perm.curve = mypval2$curve
    r1$perm.maxZs = mypval2$maxZs
    r1$perm.Z = mypval2$Z
    cat("p-value from", B, "permutations:", r1$pval.perm, "\n")
  }
  return(r1)
}

# changed interval
gseg2 = function(n, E, l0=0.05*n, l1=0.95*n, pval.appr=TRUE, skew.corr=TRUE, pval.perm=FALSE, B=100){
  l0 = ceiling(l0)
  l1 = floor(l1)
  Ebynode = vector("list", n)  
  for(i in 1:n) Ebynode[[i]]=rep(0,0)
  for(i in 1:nrow(E)){
    Ebynode[[E[i,1]]] = c(Ebynode[[E[i,1]]],E[i,2])
    Ebynode[[E[i,2]]] = c(Ebynode[[E[i,2]]],E[i,1])
  }
  temp = gcp2bynode(n,Ebynode,l0,l1)
  r1 = list(tauhat=temp$tauhat, Zmax=temp$Zmax, Z=temp$Z, R=temp$R)
  cat("Estimated change-point location:", r1$tauhat, "\n")
  cat("Zmax:", r1$Zmax, "\n")
  if (pval.appr==TRUE){
    mypval1 = pval2(n,E,Ebynode,r1$Zmax,skew.corr,l0,l1)
    r1$pval.appr = min(mypval1,1)
    cat("Approximated p-value:", r1$pval.appr, "\n")
  }
  if (pval.perm==TRUE){
    mypval2 = permpval2(n,Ebynode,r1$Zmax,B,l0,l1)
    r1$pval.perm = min(mypval2$pval,1)
    r1$perm.curve = mypval2$curve
    r1$perm.maxZs = mypval2$maxZs
    r1$perm.Z = mypval2$Z
    cat("p-value from", B, "permutations:", r1$pval.perm, "\n")
  }
  return(r1)
}


# the Nu function
Nu = function(x){ 
  y = x/2
  (1/y)*(pnorm(y)-0.5)/(y*pnorm(y) + dnorm(y))
}

# single change-point
gcp1bynode = function(n, Ebynode, n0=ceiling(0.05*n), n1=floor(0.95*n)){
# "n" is the total number of nodes.
# Ebynode[[i]] is the list of nodes that are connect to i by an edge.
# The nodes are numbered by their order in the sequence.  
# To estimate the change-point, we find the maximum of Z(t), the standardized
# version of R(t), between n1 and n2. 
  g = rep(1,n)
  R = rep(0,n)
  for(i in 1:n){
    g[i] = 0  # update g
    links = Ebynode[[i]]
    if(i==1){
      if(length(links)>0){
        R[i] = sum(rep(g[i],length(links)) != g[links])
      } else {
        R[i] = 0
      }
    } else {
      if(length(links)>0){
        add = sum(rep(g[i],length(links)) != g[links])
        subtract = length(links)-add
        R[i] = R[i-1]+add-subtract
      } else {
        R[i] = R[i-1]
      }
    }
  }
  tt = 1:n
  nodedeg = rep(0,n)
  for(i in 1:n) nodedeg[i] = length(Ebynode[[i]])
  sumEisq = sum(nodedeg^2)
  nE = sum(nodedeg)/2
  mu.t = nE* 2*tt*(n-tt)/(n*(n-1))
  p1.tt = 2*tt*(n-tt)/(n*(n-1))
  p2.tt = tt*(n-tt)*(n-2)/(n*(n-1)*(n-2))
  p3.tt = 4*tt*(n-tt)*(tt-1)*(n-tt-1)/(n*(n-1)*(n-2)*(n-3))
  A.tt = (p1.tt-2*p2.tt+p3.tt)*nE+(p2.tt-p3.tt)*sumEisq+p3.tt*nE^2
  
  Z = (mu.t-R)/sqrt(A.tt-mu.t^2)
  Z[n] = 0
  
  temp=n0:n1
  tauhat = temp[which.max(Z[n0:n1])]
  
  return(list(tauhat=tauhat,Zmax=Z[tauhat],Z=Z,R=R))
}


# changed interval
gcp2bynode = function(n, Ebynode, l0=ceiling(0.05*n), l1=floor(0.95*n)){
  Rtmp = matrix(0,n,n)
  for (i in 1:(n-1)){
    g = rep(0,n)
    for (j in (i+1):n){
      g[j] = 1 # update g
      links = Ebynode[[j]]
      if (j == (i+1)){
        if (length(links)>0){
          Rtmp[i,j] = sum(rep(g[j],length(links)) != g[links])
        }
      }else{
        if (length(links)>0){
          add = sum(rep(g[j],length(links)) != g[links])
          subtract = length(links)-add
          Rtmp[i,j] = Rtmp[i,j-1]+add-subtract
        }else{
          Rtmp[i,j] = Rtmp[i,j-1]
        }
      }
    }
  }
  tt = 1:n
  nodedeg = rep(0,n)
  for(i in 1:n) nodedeg[i] = length(Ebynode[[i]])
  sumEisq = sum(nodedeg^2)
  nE = sum(nodedeg)/2
  mu.t = nE* 2*tt*(n-tt)/(n*(n-1))
  p1.tt = 2*tt*(n-tt)/(n*(n-1))
  p2.tt = 4*tt*(n-tt)*(tt-1)*(n-tt-1)/(n*(n-1)*(n-2)*(n-3))
  V.tt = p2.tt*nE+(p1.tt/2-p2.tt)*sumEisq+(p2.tt-p1.tt^2)*nE^2

  Rv = as.vector(t(Rtmp))
  dif = matrix(0,n,n)
  for (i in 1:n){
    for (j in 1:n){
      dif[i,j] = j-i
    }
  }
  difv = as.vector(t(dif))
  ids = which(difv>0)
  Zv = rep(0,n*n)
  Zv[ids] = (mu.t[difv[ids]]-Rv[ids])/sqrt(V.tt[difv[ids]])
  ids2 = which((difv>=l0) & (difv<=l1))
  Zmax = max(Zv[ids2])
  tauhat0 = which(Zv == Zmax)
  tauhat = c(floor(tauhat0/n)+1, (tauhat0-1)%%n+1)

  return(list(tauhat=tauhat, Zmax=Zmax, Z=matrix(Zv,n,byrow=TRUE), R=Rtmp, Zv=Zv))
}

# rho_one = n h_G
rho_one = function(n, s, sumE, sumEisq){
  f1 = 4*(n-1)*(2*s*(n-s)-n)
  f2 = ((n+1)*(n-2*s)^2-2*n*(n-1))
  f3 = 4*((n-2*s)^2-n)
  f4 = 4*n*(s-1)*(n-1)*(n-s-1)
  f5 = n*(n-1)*((n-2*s)^2-(n-2))
  f6 = 4*((n-2)*(n-2*s)^2-2*s*(n-s)+n)
  n*(n-1)*(f1*sumE + f2*sumEisq - f3*sumE^2)/(2*s*(n-s)*(f4*sumE + f5*sumEisq - f6*sumE^2))
}


# p value approximation for single change-point
pval1 = function(n, E, Ebynode, Zmax, skew.corr=TRUE, n0=ceiling(0.05*n), n1=floor(0.95*n)){
  b = Zmax
  deg = rep(0,n)
  for(i in 1:n) deg[i] = length(Ebynode[[i]])
  sumE = sum(deg)/2
  sumEisq = sum(deg^2) 
  integrand = function(s){
    x = rho_one(n,s,sumE,sumEisq)
    x*Nu(sqrt(2*b^2*x))
  }
  mypval = dnorm(b)*b*integrate(integrand, n0, n1, subdivisions=3000, stop.on.error=FALSE)$value
  if (skew.corr==FALSE){
    return(mypval)
  }
  x1 = sum(deg*(deg-1))
  x2 = sum(deg*(deg-1)*(deg-2))
  x3 = 0
  for (i in 1:nrow(E)){
    x3 = x3 + (deg[E[i,1]]-1)*(deg[E[i,2]]-1)
  }  
  x4 = sum(deg*(deg-1)*(sumE-deg))
  x5 = 0
  for (i in 1:nrow(E)){
    j = E[i,1]
    k = E[i,2]
    x5 = x5 + length(which(!is.na(match(Ebynode[[j]], Ebynode[[k]]))))
  }
  s = 1:n
  x = rho_one(n,s,sumE,sumEisq)
  p1 = 2*s*(n-s)/(n*(n-1))
  p2 = 4*s*(s-1)*(n-s)*(n-s-1)/(n*(n-1)*(n-2)*(n-3))
  p3 = s*(n-s)*((n-s-1)*(n-s-2) + (s-1)*(s-2))/(n*(n-1)*(n-2)*(n-3))
  p4 = 8*s*(s-1)*(s-2)*(n-s)*(n-s-1)*(n-s-2)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
  mu = p1*sumE
  sig = sqrt(apply(cbind(p2*sumE + (p1/2-p2)*sumEisq + (p2-p1^2)*sumE^2, rep(0,n)), 1, max))  # sigma
  ER3 = p1*sumE + p1/2*3*x1 + p2*(3*sumE*(sumE-1)-3*x1) + p3*x2 + p2/2*(3*x4-6*x3) + p4*(sumE*(sumE-1)*(sumE-2)-x2-3*x4+6*x3)- 2*p4*x5
  r = (mu^3 + 3*mu*sig^2 - ER3)/sig^3
  theta_b = rep(0,n)
  pos = which(1+2*r*b>0)
  theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)
  a = x*Nu(sqrt(2*b^2*x)) * ratio
  nn = n-length(pos)
  if (nn>0.75*n){
    return(0)
  }
  if (nn>=2*(n0-1)){
    neg = which(1+2*r*b<=0)
    dif = neg[2:nn]-neg[1:(nn-1)]
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    inc = (a[id3]-a[id2])/(id3-id2)
    a[id2:1] = a[id2+1]-inc*(1:id2)
    a[(n/2+1):n] = a[(n/2):1]
    neg2 = which(a<0)
    a[neg2] = 0
  }
  integrand = function(s){
    a[s]
  }
  result = try(dnorm(b)*b*integrate(integrand, n0, n1, subdivisions=3000, stop.on.error=FALSE)$value, silent=T)
  if (is.numeric(result)){
    return(result)
  }else{
    cat("p value approximation without skewness correction is calculated.\n")
    return(mypval)
  }
}

# p value approximation for changed interval
pval2 = function(n, E, Ebynode, Zmax, skew.corr=TRUE, l0=ceiling(0.05*n), l1=floor(0.95*n)){
  b = Zmax
  deg = rep(0,n)
  for(i in 1:n) deg[i] = length(Ebynode[[i]])
  sumE = sum(deg)/2
  sumEisq = sum(deg^2) 
  integrand = function(s){
    x = rho_one(n,s,sumE,sumEisq)
    (b^2*x*Nu(sqrt(2*b^2*x)))^2*(n-s)
  }
  mypval = dnorm(b)/b*integrate(integrand, l0, l1, subdivisions=3000, stop.on.error=FALSE)$value
  if (skew.corr==FALSE){
    return(mypval)
  }
  x1 = sum(deg*(deg-1))
  x2 = sum(deg*(deg-1)*(deg-2))
  x3 = 0
  for (i in 1:nrow(E)){
    x3 = x3 + (deg[E[i,1]]-1)*(deg[E[i,2]]-1)
  }  
  x4 = sum(deg*(deg-1)*(sumE-deg))
  x5 = 0
  for (i in 1:nrow(E)){
    j = E[i,1]
    k = E[i,2]
    x5 = x5 + length(which(!is.na(match(Ebynode[[j]], Ebynode[[k]]))))
  }
  s = 1:n
  x = rho_one(n,s,sumE,sumEisq)
  p1 = 2*s*(n-s)/(n*(n-1))
  p2 = 4*s*(s-1)*(n-s)*(n-s-1)/(n*(n-1)*(n-2)*(n-3))
  p3 = s*(n-s)*((n-s-1)*(n-s-2) + (s-1)*(s-2))/(n*(n-1)*(n-2)*(n-3))
  p4 = 8*s*(s-1)*(s-2)*(n-s)*(n-s-1)*(n-s-2)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
  mu = p1*sumE
  sig = sqrt(p2*sumE + (p1/2-p2)*sumEisq + (p2-p1^2)*sumE^2)  # sigma
  ER3 = p1*sumE + p1/2*3*x1 + p2*(3*sumE*(sumE-1)-3*x1) + p3*x2 + p2/2*(3*x4-6*x3) + p4*(sumE*(sumE-1)*(sumE-2)-x2-3*x4+6*x3) - 2*p4*x5
  r = (mu^3 + 3*mu*sig^2 - ER3)/sig^3
  theta_b = rep(0,n)
  pos = which(1+2*r*b>0)
  theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)
  a = (b^2*x*Nu(sqrt(2*b^2*x)))^2 * ratio
  nn = n-length(pos)
  if (nn>0.75*n){
    return(0)
  }
  if (nn>=2*(l0-1)){
    neg = which(1+2*r*b<=0)
    dif = neg[2:nn]-neg[1:(nn-1)]
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    inc = (a[id3]-a[id2])/ceiling(0.03*n)
    a[id2:1] = a[id2+1]-inc*(1:id2)
    a[(n/2+1):n] = a[(n/2):1]
    neg2 = which(a<0)
    a[neg2] = 0
  }
  integrand = function(s){
    a[s]*(n-s)
  }
  result = try(dnorm(b)/b*integrate(integrand, l0, l1, subdivisions=3000, stop.on.error=FALSE)$value, silent=T)
  if (is.numeric(result)){
    return(result)
  }else{
    cat("p value approximation without skewness correction is reported.\n")
    return(mypval)
  }
}

# p value from permutation for single change point
permpval1 = function(n, Ebynode, Zmax, B=100, n0=ceiling(0.05*n), n1=floor(0.95*n)){
# Computes the pvalue P(max_{n1<=t<=n2} Z(t)>b) by permuting the nodes in the graph.
   Z = matrix(0,B,n)
  for(b in 1:B){
    if(b%%1000 ==0) {
      cat(b, "permutations completed.\n")
    }
    perm = sample(n)
    permmatch = rep(0,n)
    for(i in 1:n) permmatch[perm[i]] = i
    Ebnstar =  vector("list", n)
    for(i in 1:n){
      oldlinks = Ebynode[[permmatch[i]]]
      Ebnstar[[i]] = perm[oldlinks]
    }
    gcpstar=gcp1bynode(n,Ebnstar,n0,n1)
    Z[b,] = gcpstar$Z
  }
  
  maxZ = apply(Z[,n0:n1],1,max)
  maxZs=sort(maxZ)
  p=1-(0:(B-1))/B
  ## if (isplot){
  ##   plot(maxZs,p,type="l")
  ## }
  return(list(pval=length(which(maxZs>=Zmax))/B, curve=cbind(maxZs,p), maxZs=maxZs, Z=Z))
}

# p value from permutation for changed interval
permpval2 = function(n,Ebynode,Zmax,B=100,l0=ceiling(0.05*n),l1=floor(0.95*n)){
# Computes the pvalue for changed interval by permuting the nodes in the graph.  
  Z = matrix(nrow=B,ncol=n*n)
  for(b in 1:B){
    if(b%%100 ==0) {
      cat(b, "permutations completed.\n")
    }
    perm = sample(n)
    permmatch = rep(0,n)
    for(i in 1:n) permmatch[perm[i]] = i
    Ebnstar =  vector("list", n)
    for(i in 1:n){
      oldlinks = Ebynode[[permmatch[i]]]
      Ebnstar[[i]] = perm[oldlinks]
    }
    gcpstar=gcp2bynode(n,Ebnstar,l0,l1)
    Z[b,] = gcpstar$Zv
  }
  dif = matrix(0,n,n)
  for (i in 1:n){
    for (j in 1:n){
      dif[i,j] = j-i
    }
  }
  difv = as.vector(t(dif))
  ids = which((difv>=l0)+(difv<=l1)==2)
  maxZ = apply(Z[,ids],1,max)
  maxZs = sort(maxZ)
  p=1-(0:(B-1))/B
  ## if (isplot){
  ##   plot(maxZs,p,type="l")
  ## }
  return(list(pval=length(which(maxZs>=Zmax))/B, curve=cbind(maxZs,p), maxZs=maxZs, Z=Z))
}




