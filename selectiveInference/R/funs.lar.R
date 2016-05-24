# We compute the least angle regression (LAR) path given
# a response vector y and predictor matrix x.  We assume
# that x has columns in general position.

# NOTE: the df estimates at each lambda_k can be thought of as the df
# for all solutions corresponding to lambda in (lambda_k,lambda_{k-1}),
# the open interval to the *right* of the current lambda_k.

# NOTE: x having columns in general position implies that the
# centered x satisfies a modified version of the general position
# condition, where we replace k < min(n,p) by k < min(n-1,p) in
# the definition. This is still sufficient to imply the uniqueness
# of the lasso solution, on the centered x

lar <- function(x, y, maxsteps=2000, minlam=0, intercept=TRUE, normalize=TRUE,
                verbose=FALSE) {

  this.call = match.call()
  checkargs.xy(x=x,y=y)
  
  # Center and scale, etc.
  obj = standardize(x,y,intercept,normalize)
  x = obj$x
  y = obj$y
  bx = obj$bx
  by = obj$by
  sx = obj$sx
  n = nrow(x)
  p = ncol(x)

  #####
  # Find the first variable to enter and its sign
  uhat = t(x)%*%y
  ihit = which.max(abs(uhat))   # Hitting coordinate
  hit = abs(uhat[ihit])         # Critical lambda
  s = Sign(uhat[ihit])          # Sign

  if (verbose) {
    cat(sprintf("1. lambda=%.3f, adding variable %i, |A|=%i...",
                hit,ihit,1))
  }

  # Now iteratively find the new LAR estimate, and
  # the next critical lambda

  # Things to keep track of, and return at the end
  buf = min(maxsteps,500)
  lambda = numeric(buf)      # Critical lambdas
  action = numeric(buf)      # Actions taken
  df = numeric(buf)          # Degrees of freedom
  beta = matrix(0,p,buf)     # LAR estimates
  
  lambda[1] = hit
  action[1] = ihit
  df[1] = 0
  beta[,1] = 0

  # Gamma matrix!
  gbuf = max(2*p*3,2000)     # Space for 3 steps, at least
  gi = 0
  Gamma = matrix(0,gbuf,n)
  Gamma[gi+Seq(1,p-1),] = t(s*x[,ihit]+x[,-ihit]); gi = gi+p-1
  Gamma[gi+Seq(1,p-1),] = t(s*x[,ihit]-x[,-ihit]); gi = gi+p-1
  Gamma[gi+1,] = t(s*x[,ihit]); gi = gi+1

  # nk, regression contrast, M plus
  nk = mp = numeric(buf)
  vreg = matrix(0,buf,n)

  nk[1] = gi
  vreg[1,] = s*x[,ihit] / sum(x[,ihit]^2)
  if (p > 1) {
    c = t(as.numeric(Sign(t(x)%*%y)) * t(x))
    ratio = t(c[,-ihit])%*%c[,ihit]/sum(c[,ihit]^2)
    ip = 1-ratio > 0
    crit = (t(c[,-ihit])%*%y - ratio*sum(c[,ihit]*y))/(1-ratio)
    mp[1] = max(max(crit[ip]),0)
  }

  # Other things to keep track of, but not return
  r = 1                      # Size of active set
  A = ihit                   # Active set
  I = Seq(1,p)[-ihit]        # Inactive set
  X1 = x[,ihit,drop=FALSE]   # Matrix X[,A]
  X2 = x[,-ihit,drop=FALSE]  # Matrix X[,I]
  k = 2                      # What step are we at?

  # Compute a skinny QR decomposition of X1
  obj = qr(X1)
  Q = qr.Q(obj,complete=TRUE)
  Q1 = Q[,1,drop=FALSE];
  Q2 = Q[,-1,drop=FALSE]
  R = qr.R(obj)
  
  # Throughout the algorithm, we will maintain
  # the decomposition X1 = Q1*R. Dimenisons:
  # X1: n x r
  # Q1: n x r
  # Q2: n x (n-r)
  # R:  r x r
    
  while (k<=maxsteps && lambda[k-1]>=minlam) {
    ##########
    # Check if we've reached the end of the buffer
    if (k > length(lambda)) {
      buf = length(lambda)
      lambda = c(lambda,numeric(buf))
      action = c(action,numeric(buf))
      df = c(df,numeric(buf))
      beta = cbind(beta,matrix(0,p,buf))
      nk = c(nk,numeric(buf))
      mp = c(mp,numeric(buf))
      vreg = rbind(vreg,matrix(0,buf,n))
    }

    # Key quantities for the hitting times
    a = backsolve(R,t(Q1)%*%y)
    b = backsolve(R,backsolve(R,s,transpose=TRUE))
    aa = as.numeric(t(X2) %*% (y - X1 %*% a))
    bb = as.numeric(t(X2) %*% (X1 %*% b))
    
    # If the inactive set is empty, nothing will hit
    if (r==min(n-intercept,p)) hit = 0

    # Otherwise find the next hitting time
    else {
      shits = Sign(aa)
      hits = aa/(shits-bb)

      # Make sure none of the hitting times are larger
      # than the current lambda 
      hits[hits>lambda[k-1]] = 0
        
      ihit = which.max(hits)
      hit = hits[ihit]
      shit = shits[ihit]
    }

    # Stop if the next critical point is negative
    if (hit<=0) break
    
    # Record the critical lambda and solution
    lambda[k] = hit
    action[k] = I[ihit]
    df[k] = r
    beta[A,k] = a-hit*b
        
    # Gamma matrix!
    if (gi + 2*p > nrow(Gamma)) Gamma = rbind(Gamma,matrix(0,2*p+gbuf,n))
    X2perp = X2 - X1 %*% backsolve(R,t(Q1)%*%X2)
    c = t(t(X2perp)/(shits-bb))
    Gamma[gi+Seq(1,p-r),] = shits*t(X2perp); gi = gi+p-r
    Gamma[gi+Seq(1,p-r-1),] = t(c[,ihit]-c[,-ihit]); gi = gi+p-r-1
    Gamma[gi+1,] = t(c[,ihit]); gi = gi+1

    # nk, regression contrast, M plus
    nk[k] = gi
    vreg[k,] = shit*X2perp[,ihit] / sum(X2perp[,ihit]^2)
    if (ncol(c) > 1) {
      ratio = t(c[,-ihit])%*%c[,ihit]/sum(c[,ihit]^2)
      ip = 1-ratio > 0
      crit = (t(c[,-ihit])%*%y - ratio*sum(c[,ihit]*y))/(1-ratio)
      mp[k] = max(max(crit[ip]),0)
    }
    
    # Update all of the variables
    r = r+1
    A = c(A,I[ihit])
    I = I[-ihit]
    s = c(s,shit)
    X1 = cbind(X1,X2[,ihit])
    X2 = X2[,-ihit,drop=FALSE]

    # Update the QR decomposition
    obj = updateQR(Q1,Q2,R,X1[,r])
    Q1 = obj$Q1
    Q2 = obj$Q2
    R = obj$R
     
    if (verbose) {
      cat(sprintf("\n%i. lambda=%.3f, adding variable %i, |A|=%i...",
                  k,hit,A[r],r))
    }
            
    # Step counter
    k = k+1
  }

  # Trim
  lambda = lambda[Seq(1,k-1)]
  action = action[Seq(1,k-1)]
  df = df[Seq(1,k-1),drop=FALSE]
  beta = beta[,Seq(1,k-1),drop=FALSE]
  Gamma = Gamma[Seq(1,gi),,drop=FALSE]
  nk = nk[Seq(1,k-1)]
  mp = mp[Seq(1,k-1)]
  vreg = vreg[Seq(1,k-1),,drop=FALSE]
  
  # If we reached the maximum number of steps
  if (k>maxsteps) {
    if (verbose) {
      cat(sprintf("\nReached the maximum number of steps (%i),",maxsteps))
      cat(" skipping the rest of the path.")
    }
    completepath = FALSE
    bls = NULL
  }

  # If we reached the minimum lambda
  else if (lambda[k-1]<minlam) {
    if (verbose) {
      cat(sprintf("\nReached the minimum lambda (%.3f),",minlam))
      cat(" skipping the rest of the path.")
    }
    completepath = FALSE
    bls = NULL
  }
  
  # Otherwise, note that we completed the path
  else {
    completepath = TRUE
    
    # Record the least squares solution. Note that
    # we have already computed this
    bls = rep(0,p)
    bls[A] = a
  }

  if (verbose) cat("\n")
  
  # Adjust for the effect of centering and scaling
  if (intercept) df = df+1
  if (normalize) beta = beta/sx
  if (normalize && completepath) bls = bls/sx
  
  # Assign column names
  colnames(beta) = as.character(round(lambda,3))

  out = list(lambda=lambda,action=action,sign=s,df=df,beta=beta,
    completepath=completepath,bls=bls,
    Gamma=Gamma,nk=nk,vreg=vreg,mp=mp,x=x,y=y,bx=bx,by=by,sx=sx,
    intercept=intercept,normalize=normalize,call=this.call) 
  class(out) = "lar"
  return(out)
}

##############################

# Downdate the QR factorization, after a column has
# been deleted. Here Q1 is m x n, Q2 is m x k, and
# R is n x n.

downdateQR <- function(Q1,Q2,R,col) {
  m = nrow(Q1)
  n = ncol(Q1)
  
  a = .C("downdate1",
    Q1=as.double(Q1),
    R=as.double(R),
    col=as.integer(col-1),
    m=as.integer(m),
    n=as.integer(n),
    dup=FALSE,
    package="selectiveInference")

  Q1 = matrix(a$Q1,nrow=m)
  R = matrix(a$R,nrow=n)

  # Re-structure: add a column to Q2, delete one from
  # Q1, and trim R
  Q2 = cbind(Q2,Q1[,n])
  Q1 = Q1[,-n,drop=FALSE]
  R = R[-n,-col,drop=FALSE]

  return(list(Q1=Q1,Q2=Q2,R=R))
}

# Update the QR factorization, after a column has been
# added. Here Q1 is m x n, Q2 is m x k, and R is n x n.

updateQR <- function(Q1,Q2,R,col) {
  m = nrow(Q1)
  n = ncol(Q1)
  k = ncol(Q2)
  
  a = .C("update1",
    Q2=as.double(Q2),
    w=as.double(t(Q2)%*%col),
    m=as.integer(m),
    k=as.integer(k),
    dup=FALSE,
    package="selectiveInference")

  Q2 = matrix(a$Q2,nrow=m)
  w = c(t(Q1)%*%col,a$w)

  # Re-structure: delete a column from Q2, add one to
  # Q1, and expand R
  Q1 = cbind(Q1,Q2[,1])
  Q2 = Q2[,-1,drop=FALSE]
  R = rbind(R,rep(0,n))
  R = cbind(R,w[Seq(1,n+1)])

  return(list(Q1=Q1,Q2=Q2,R=R))
}

##############################

# Coefficient function for lar

coef.lar <- function(object, s, mode=c("step","lambda"), ...) {
  mode = match.arg(mode)

  if (object$completepath) {
    k = length(object$action)+1
    lambda = c(object$lambda,0)
    beta = cbind(object$beta,object$bls)
  } else {
    k = length(object$action)
    lambda = object$lambda
    beta = object$beta
  }
  
  if (mode=="step") {
    if (min(s)<0 || max(s)>k) stop(sprintf("s must be between 0 and %i",k))
    knots = 1:k
    dec = FALSE
  } else {
    if (min(s)<min(lambda)) stop(sprintf("s must be >= %0.3f",min(lambda)))
    knots = lambda
    dec = TRUE
  }
  
  return(coef.interpolate(beta,s,knots,dec))
}

# Prediction function for lar

predict.lar <- function(object, newx, s, mode=c("step","lambda"), ...) {
  beta = coef.lar(object,s,mode)
  if (missing(newx)) newx = scale(object$x,FALSE,1/object$sx)
  else newx = scale(newx,object$bx,FALSE)
  return(newx %*% beta + object$by)
}

coef.lasso <- coef.lar
predict.lasso <- predict.lar

##############################

# Lar inference function

larInf <- function(obj, sigma=NULL, alpha=0.1, k=NULL, type=c("active","all","aic"), 
                   gridrange=c(-100,100), bits=NULL, mult=2, ntimes=2, verbose=FALSE) {
  
  this.call = match.call()
  type = match.arg(type)
  checkargs.misc(sigma=sigma,alpha=alpha,k=k,
                 gridrange=gridrange,mult=mult,ntimes=ntimes)
  if (class(obj) != "lar") stop("obj must be an object of class lar")
  if (is.null(k) && type=="active") k = length(obj$action)
  if (is.null(k) && type=="all") stop("k must be specified when type = all")
  if (!is.null(bits) && !requireNamespace("Rmpfr",quietly=TRUE)) {
    warning("Package Rmpfr is not installed, reverting to standard precision")
    bits = NULL
  }
  
  k = min(k,length(obj$action)) # Round to last step
  x = obj$x
  y = obj$y
  p = ncol(x)
  n = nrow(x)
  G = obj$Gamma
  nk = obj$nk
  sx = obj$sx

  if (is.null(sigma)) {
    if (n >= 2*p) {
      oo = obj$intercept
      sigma = sqrt(sum(lsfit(x,y,intercept=oo)$res^2)/(n-p-oo))
    }
    else {
      sigma = sd(y)
      warning(paste(sprintf("p > n/2, and sd(y) = %0.3f used as an estimate of sigma;",sigma),
                    "you may want to use the estimateSigma function"))
    }
  }
  
  pv.spacing = pv.modspac = pv.covtest = khat = NULL
  
  if (type == "active") {
    pv = vlo = vup = numeric(k) 
    vmat = matrix(0,k,n)
    ci = tailarea = matrix(0,k,2)
    pv.spacing = pv.modspac = pv.covtest = numeric(k)
    vreg = obj$vreg[1:k,,drop=FALSE]
    sign = obj$sign[1:k]
    vars = obj$action[1:k]

    for (j in 1:k) {
      if (verbose) cat(sprintf("Inference for variable %i ...\n",vars[j]))
      
      Gj = G[1:nk[j],]
      uj = rep(0,nk[j])
      vj = vreg[j,]
      mj = sqrt(sum(vj^2))
      vj = vj / mj              # Standardize (divide by norm of vj)
      a = poly.pval(y,Gj,uj,vj,sigma,bits)
      pv[j] = a$pv
      sxj = sx[vars[j]]
      vlo[j] = a$vlo * mj / sxj # Unstandardize (mult by norm of vj / sxj)
      vup[j] = a$vup * mj / sxj # Unstandardize (mult by norm of vj)
      vmat[j,] = vj * mj / sxj  # Unstandardize (mult by norm of vj / sxj)
    
      a = poly.int(y,Gj,uj,vj,sigma,alpha,gridrange=gridrange,
        flip=(sign[j]==-1),bits=bits)
      ci[j,] = a$int * mj / sxj # Unstandardize (mult by norm of vj / sxj) 
      tailarea[j,] = a$tailarea
      
      pv.spacing[j] = spacing.pval(obj,sigma,j)
      pv.modspac[j] = modspac.pval(obj,sigma,j)
      pv.covtest[j] = covtest.pval(obj,sigma,j)
    }

    khat = forwardStop(pv,alpha)
  }
  
  else {
    if (type == "aic") {
      out = aicStop(x,y,obj$action[1:k],obj$df[1:k],sigma,mult,ntimes)
      khat = out$khat
      m = out$stopped * ntimes
      G = rbind(out$G,G[1:nk[khat+m],])  # Take ntimes more steps past khat
      u = c(out$u,rep(0,nk[khat+m]))     # (if we need to)
      kk = khat
    }
    else {
      G = G[1:nk[k],]
      u = rep(0,nk[k])
      kk = k
    }
    
    pv = vlo = vup = numeric(kk) 
    vmat = matrix(0,kk,n)
    ci = tailarea = matrix(0,kk,2)
    sign = numeric(kk)
    vars = obj$action[1:kk]
    xa = x[,vars]
    M = pinv(crossprod(xa)) %*% t(xa)
    
    for (j in 1:kk) {
      if (verbose) cat(sprintf("Inference for variable %i ...\n",vars[j]))
      
      vj = M[j,]
      mj = sqrt(sum(vj^2))
      vj = vj / mj             # Standardize (divide by norm of vj)
      sign[j] = sign(sum(vj*y)) 
      vj = sign[j] * vj
      Gj = rbind(G,vj)
      uj = c(u,0)

      a = poly.pval(y,Gj,uj,vj,sigma,bits)
      pv[j] = a$pv
      sxj = sx[vars[j]]
      vlo[j] = a$vlo * mj / sxj # Unstandardize (mult by norm of vj / sxj)
      vup[j] = a$vup * mj / sxj # Unstandardize (mult by norm of vj / sxj)
      vmat[j,] = vj * mj / sxj  # Unstandardize (mult by norm of vj / sxj)

      a = poly.int(y,Gj,uj,vj,sigma,alpha,gridrange=gridrange,
        flip=(sign[j]==-1),bits=bits)
      ci[j,] = a$int * mj / sxj # Unstandardize (mult by norm of vj / sxj)
      tailarea[j,] = a$tailarea
    }
  }
  
  out = list(type=type,k=k,khat=khat,pv=pv,ci=ci,
    tailarea=tailarea,vlo=vlo,vup=vup,vmat=vmat,y=y,
    pv.spacing=pv.spacing,pv.modspac=pv.modspac,pv.covtest=pv.covtest,
    vars=vars,sign=sign,sigma=sigma,alpha=alpha,
    call=this.call)
  class(out) = "larInf"
  return(out)
}

##############################

spacing.pval <- function(obj, sigma, k) {
  v = obj$Gamma[obj$nk[k],]
  sd = sigma*sqrt(sum(v^2))
  a = obj$mp[k]
  
  if (k==1) b = Inf
  else b = obj$lambda[k-1]
  
  return(tnorm.surv(obj$lambda[k],0,sd,a,b))
}

modspac.pval <- function(obj, sigma, k) {
  v = obj$Gamma[obj$nk[k],]
  sd = sigma*sqrt(sum(v^2))

  if (k < length(obj$action)) a = obj$lambda[k+1]
  else if (obj$completepath) a = 0
  else {
    warning(sprintf("Modified spacing p-values at step %i require %i steps of the lar path",k,k+1))
    return(NA)
  }
      
  if (k==1) b = Inf
  else b = obj$lambda[k-1]

  return(tnorm.surv(obj$lambda[k],0,sd,a,b))
}

covtest.pval <- function(obj, sigma, k) {
  A = which(obj$beta[,k]!=0)
  sA = sign(obj$beta[A,k])
  lam1 = obj$lambda[k]
  j = obj$action[k]

  if (k < length(obj$action)) {
    lam2 = obj$lambda[k+1]
    sj = sign(obj$beta[j,k+1])
  } else if (obj$completepath) {
    lam2 = 0
    sj = sign(obj$bls[j])
  } else {
    warning(sprintf("Cov test p-values at step %i require %i steps of the lar path",k,k+1))
    return(NA)
  }

  x = obj$x
  if (length(A)==0) term1 = 0
  else term1 = x[,A,drop=F] %*% solve(crossprod(x[,A,drop=F]),sA)
  term2 = x[,c(A,j),drop=F] %*% solve(crossprod(x[,c(A,j),drop=F]),c(sA,sj))
  c = sum((term2 - term1)^2)
  t = c * lam1 * (lam1-lam2) / sigma^2
  return(1-pexp(t))
}

##############################

print.lar <- function(x, ...) {
  cat("\nCall:\n")
  dput(x$call)
  
  cat("\nSequence of LAR moves:\n")
  nsteps = length(x$action)
  tab = cbind(1:nsteps,x$action,x$sign)
  colnames(tab) = c("Step","Var","Sign")
  rownames(tab) = rep("",nrow(tab))
  print(tab)
  invisible()
}

print.larInf <- function(x, tailarea=TRUE, ...) {
  cat("\nCall:\n")
  dput(x$call)

  cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n",
              x$sigma))

  if (x$type == "active") {
    cat(sprintf("\nSequential testing results with alpha = %0.3f\n",x$alpha))
    cat("",fill=T)
    tab = cbind(1:length(x$pv),x$vars,
      round(x$sign*x$vmat%*%x$y,3),
      round(x$sign*x$vmat%*%x$y/(x$sigma*sqrt(rowSums(x$vmat^2))),3),
      round(x$pv,3),round(x$ci,3),round(x$pv.spacing,3),round(x$pv.cov,3)) 
    colnames(tab) = c("Step", "Var", "Coef", "Z-score", "P-value",
              "LowConfPt", "UpConfPt", "Spacing", "CovTest")
    if (tailarea) {
      tab = cbind(tab,round(x$tailarea,3))
      colnames(tab)[(ncol(tab)-1):ncol(tab)] = c("LowTailArea","UpTailArea")
    }
    rownames(tab) = rep("",nrow(tab))
    print(tab)

    cat(sprintf("\nEstimated stopping point from ForwardStop rule = %i\n",x$khat))
  }

  else if (x$type == "all") {
    cat(sprintf("\nTesting results at step = %i, with alpha = %0.3f\n",x$k,x$alpha))
    cat("",fill=T)
    tab = cbind(x$vars,
      round(x$sign*x$vmat%*%x$y,3),
      round(x$sign*x$vmat%*%x$y/(x$sigma*sqrt(rowSums(x$vmat^2))),3),
      round(x$pv,3),round(x$ci,3))
    colnames(tab) = c("Var", "Coef", "Z-score", "P-value", "LowConfPt", "UpConfPt")
    if (tailarea) {
      tab = cbind(tab,round(x$tailarea,3))
      colnames(tab)[(ncol(tab)-1):ncol(tab)] = c("LowTailArea","UpTailArea")
    }
    rownames(tab) = rep("",nrow(tab))
    print(tab)
  }

  else if (x$type == "aic") {
    cat(sprintf("\nTesting results at step = %i, with alpha = %0.3f\n",x$khat,x$alpha))
     cat("",fill=T)
    tab = cbind(x$vars,
      round(x$sign*x$vmat%*%x$y,3),
      round(x$sign*x$vmat%*%x$y/(x$sigma*sqrt(rowSums(x$vmat^2))),3),
      round(x$pv,3),round(x$ci,3))
    colnames(tab) = c("Var", "Coef", "Z-score", "P-value", "LowConfPt", "UpConfPt")
    if (tailarea) {
      tab = cbind(tab,round(x$tailarea,3))
      colnames(tab)[(ncol(tab)-1):ncol(tab)] = c("LowTailArea","UpTailArea")
    }
    rownames(tab) = rep("",nrow(tab))
    print(tab)
    
    cat(sprintf("\nEstimated stopping point from AIC rule = %i\n",x$khat))
  }

  invisible()
}

plot.lar <- function(x, xvar=c("norm","step","lambda"), breaks=TRUE,
                     omit.zeros=TRUE, var.labels=TRUE, ...) {
  
  if (x$completepath) {
    k = length(x$action)+1
    lambda = c(x$lambda,0)
    beta = cbind(x$beta,x$bls)
  } else {
    k = length(x$action)
    lambda = x$lambda
    beta = x$beta
  }
  p = nrow(beta)
  
  xvar = match.arg(xvar)
  if (xvar=="norm") {
    xx = colSums(abs(beta))
    xlab = "L1 norm"
  } else if (xvar=="step") {
    xx = 1:k
    xlab = "Step"
  } else {
    xx = lambda
    xlab = "Lambda"
  }

  if (omit.zeros) {
    good.inds = matrix(FALSE,p,k)
    good.inds[beta!=0] = TRUE
    changes = t(apply(beta,1,diff))!=0
    good.inds[cbind(changes,rep(F,p))] = TRUE
    good.inds[cbind(rep(F,p),changes)] = TRUE
    beta[!good.inds] = NA
  }

  plot(c(),c(),xlim=range(xx,na.rm=T),ylim=range(beta,na.rm=T),
       xlab=xlab,ylab="Coefficients",main="Least angle regression path",...)
  abline(h=0,lwd=2)
  matplot(xx,t(beta),type="l",lty=1,add=TRUE)
  if (breaks) abline(v=xx,lty=2)
  if (var.labels) axis(4,at=beta[,k],labels=1:p,cex=0.8,adj=0) 
  invisible()
}
