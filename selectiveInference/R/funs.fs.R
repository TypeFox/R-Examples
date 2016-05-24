# We compute the forward stepwise regression (FS) path given 
# a response vector y and predictor matrix x.  We assume
# that x has columns in general position.

fs <- function(x, y, maxsteps=2000, intercept=TRUE, normalize=TRUE,
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
  # To keep consistent with the lar function, we parametrize
  # so that the first step has all zero coefficients,
  # Also, an interesting note: the effective "lambda" (maximal
  # correlation with the residual) may increase with stepwise!
  # So we don't keep track of it

  #####
  # Find the first variable to enter and its sign
  working_x = scale(x,center=F,scale=sqrt(colSums(x^2)))
  score = t(working_x)%*%y
  i_hit = which.max(abs(score))   # Hitting coordinate
  sign_hit = Sign(score[i_hit])   # Sign
  signs = sign_hit                # later signs will be appended to `signs`

  if (verbose) {
    cat(sprintf("1. Adding variable %i, |A|=%i...",i_hit,1))
  }
  
  # Now iteratively find the new FS estimates

  # Things to keep track of, and return at the end
  # JT: I guess the "buf" just saves us from making huge
  # matrices we don't need?

  buf = min(maxsteps,500)
  action = numeric(buf)      # Actions taken
  df = numeric(buf)          # Degrees of freedom
  beta = matrix(0,p,buf)     # FS estimates
  
  action[1] = i_hit
  df[1] = 0
  beta[,1] = 0

  # Gamma matrix!
  gbuf = max(2*p*3,2000)     # Space for 3 steps, at least
  gi = 0                     # index into rows of Gamma matrix

  Gamma = matrix(0,gbuf,n)
  Gamma[gi+Seq(1,p-1),] = t(sign_hit*working_x[,i_hit]+working_x[,-i_hit]); gi = gi+p-1
  Gamma[gi+Seq(1,p-1),] = t(sign_hit*working_x[,i_hit]-working_x[,-i_hit]); gi = gi+p-1
  Gamma[gi+1,] = t(sign_hit*working_x[,i_hit]); gi = gi+1

  # nconstraint
  nconstraint = numeric(buf)
  vreg = matrix(0,buf,n)
  nconstraint[1] = gi
  vreg[1,] = sign_hit*x[,i_hit] / sum(x[,i_hit]^2)

  # Other things to keep track of, but not return
  r = 1                      # Size of active set
  A = i_hit                  # Active set -- JT: isn't this basically the same as action?
  I = Seq(1,p)[-i_hit]       # Inactive set
  X_active = x[,i_hit,drop=FALSE]   # Matrix X[,A]
  X_inactive = x[,-i_hit,drop=FALSE]  # Matrix X[,I]
  k = 2                      # What step are we at?
                             # JT Why keep track of r and k instead of just saying k=r+1?

  # Compute a skinny QR decomposition of X_active
  # JT: obs was used as variable name above -- this is something different, no?
  # changed it to qr_X

  qr_X = qr(X_active)
  Q = qr.Q(qr_X,complete=TRUE)
  Q_active = Q[,1,drop=FALSE];
  Q_inactive = Q[,-1,drop=FALSE]
  R = qr.R(qr_X)
  
  # Throughout the algorithm, we will maintain
  # the decomposition X_active = Q_active*R. Dimensions:
  # X_active: n x r
  # Q_active: n x r
  # Q_inactive: n x (n-r)
  # R:  r x r
    
  while (k<=maxsteps) {
    ##########
    # Check if we've reached the end of the buffer
    if (k > length(action)) {
      buf = length(action)
      action = c(action,numeric(buf))
      df = c(df,numeric(buf))
      beta = cbind(beta,matrix(0,p,buf))
      nconstraint = c(nconstraint,numeric(buf))
      vreg = rbind(vreg,matrix(0,buf,n))
    }

    # Key quantities for the next entry
    keepLs=backsolve(R,t(Q_active)%*%X_inactive)
    X_inactive_resid = X_inactive - X_active %*% keepLs
    working_x = scale(X_inactive_resid,center=F,scale=sqrt(colSums(X_inactive_resid^2)))
    score = as.numeric(t(working_x)%*%y)
    
    # If the inactive set is empty, nothing will hit
    if (r==min(n-intercept,p)) break

    # Otherwise find the next hitting time
    else {
      sign_score = Sign(score)
      abs_score = sign_score * score
      i_hit = which.max(abs_score)
      sign_hit = sign_score[i_hit]
    }
    
    # Record the solution
    # what is the difference between "action" and "A"?

    action[k] = I[i_hit] 
    df[k] = r
    beta[A,k] = backsolve(R,t(Q_active)%*%y)
        
    # Gamma matrix!
    if (gi + 2*p > nrow(Gamma)) Gamma = rbind(Gamma,matrix(0,2*p+gbuf,n))
    working_x = t(sign_score*t(working_x))
    Gamma[gi+Seq(1,p-r),] = t(working_x); gi = gi+p-r
    Gamma[gi+Seq(1,p-r-1),] = t(working_x[,i_hit]-working_x[,-i_hit]); gi = gi+p-r-1
    Gamma[gi+1,] = t(working_x[,i_hit]); gi = gi+1

    # nconstraint, regression contrast
    nconstraint[k] = gi
    vreg[k,] = sign_hit*X_inactive_resid[,i_hit] / sum(X_inactive_resid[,i_hit]^2)

    # Update all of the variables
    r = r+1
    A = c(A,I[i_hit])
    I = I[-i_hit]
    signs = c(signs,sign_hit)
    X_active = cbind(X_active,X_inactive[,i_hit])
    X_inactive = X_inactive[,-i_hit,drop=FALSE]

    # Update the QR decomposition
    updated_qr = updateQR(Q_active,Q_inactive,R,X_active[,r])
    Q_active = updated_qr$Q1

    # JT: why do we store Q_inactive? Doesn't seem to be used.
    Q_inactive = updated_qr$Q2
    R = updated_qr$R
     
    if (verbose) {
      cat(sprintf("\n%i. Adding variable %i, |A|=%i...",k,A[r],r))
    }
            
    # Step counter
    k = k+1
  }

  # Trim
  action = action[Seq(1,k-1)]
  df = df[Seq(1,k-1),drop=FALSE]
  beta = beta[,Seq(1,k-1),drop=FALSE]
  Gamma = Gamma[Seq(1,gi),,drop=FALSE]
  nconstraint = nconstraint[Seq(1,k-1)]
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

  # Otherwise, note that we completed the path
  else {
    completepath = TRUE
    
    # Record the least squares solution. Note that
    # we have already computed this
    bls = rep(0,p)
   if(length(keepLs)>0)  bls[A] = keepLs
  }

  if (verbose) cat("\n")
  
  # Adjust for the effect of centering and scaling
  if (intercept) df = df+1
  if (normalize) beta = beta/sx
  if (normalize && completepath) bls = bls/sx
  
  # Assign column names
  colnames(beta) = as.character(Seq(1,k-1))

  out = list(action=action,sign=signs,df=df,beta=beta,
    completepath=completepath,bls=bls,
    Gamma=Gamma,nconstraint=nconstraint,vreg=vreg,x=x,y=y,bx=bx,by=by,sx=sx,
    intercept=intercept,normalize=normalize,call=this.call) 
  class(out) = "fs"
  return(out)
}

##############################

# Coefficient function for fs

coef.fs <- function(object, s, ...) {
  if (object$completepath) {
    k = length(object$action)+1
    beta = cbind(object$beta,object$bls)
  } else {
    k = length(object$action)
    beta = object$beta
  }
  
  if (min(s)<0 || max(s)>k) stop(sprintf("s must be between 0 and %i",k))
  knots = 1:k
  dec = FALSE
  return(coef.interpolate(beta,s,knots,dec))
}

# Prediction function for fs

predict.fs <- function(object, newx, s, ...) {
  beta = coef.fs(object,s)
  if (missing(newx)) newx = scale(object$x,FALSE,1/object$sx)
  else newx = scale(newx,object$bx,FALSE)
  return(newx %*% beta + object$by)
}

##############################

# FS inference function

fsInf <- function(obj, sigma=NULL, alpha=0.1, k=NULL, type=c("active","all","aic"), 
                  gridrange=c(-100,100), bits=NULL, mult=2, ntimes=2, verbose=FALSE) {
  
  this.call = match.call()
  type = match.arg(type)
  checkargs.misc(sigma=sigma,alpha=alpha,k=k,
                 gridrange=gridrange,mult=mult,ntimes=ntimes)
  if (class(obj) != "fs") stop("obj must be an object of class fs")
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
  nconstraint = obj$nconstraint
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

  khat = NULL
  
  if (type == "active") {
    pv = vlo = vup = numeric(k) 
    vmat = matrix(0,k,n)
    ci = tailarea = matrix(0,k,2)
    vreg = obj$vreg[1:k,,drop=FALSE]
    sign = obj$sign[1:k]
    vars = obj$action[1:k]

    for (j in 1:k) {
      if (verbose) cat(sprintf("Inference for variable %i ...\n",vars[j]))

      Gj = G[1:nconstraint[j],]
      uj = rep(0,nconstraint[j])
      vj = vreg[j,]
      mj = sqrt(sum(vj^2)) 
      vj = vj / mj              # Standardize (divide by norm of vj)
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

    khat = forwardStop(pv,alpha)
  }
  
  else {
    if (type == "aic") {
      out = aicStop(x,y,obj$action[1:k],obj$df[1:k],sigma,mult,ntimes)
      khat = out$khat
      m = out$stopped * ntimes
      G = rbind(out$G,G[1:nconstraint[khat+m],])  # Take ntimes more steps past khat
      u = c(out$u,rep(0,nconstraint[khat+m]))     # (if we need to)
      kk = khat
    }
    else {
      G = G[1:nconstraint[k],]
      u = rep(0,nconstraint[k])
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
      vj = vj / mj              # Standardize (divide by norm of vj)
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
  
  # JT: why do we output vup, vlo? Are they used somewhere else?

  out = list(type=type,k=k,khat=khat,pv=pv,ci=ci,
    tailarea=tailarea,vlo=vlo,vup=vup,vmat=vmat,y=y,
    vars=vars,sign=sign,sigma=sigma,alpha=alpha,
    call=this.call)
  class(out) = "fsInf"
  return(out)
}

##############################


##############################

print.fs <- function(x, ...) {
  cat("\nCall:\n")
  dput(x$call)
  
  cat("\nSequence of FS moves:\n")
  nsteps = length(x$action)
  tab = cbind(1:nsteps,x$action,x$sign)
  colnames(tab) = c("Step","Var","Sign")
  rownames(tab) = rep("",nrow(tab))
  print(tab)
  invisible()
}

print.fsInf <- function(x, tailarea=TRUE, ...) {
  cat("\nCall:\n")
  dput(x$call)

  cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n",
              x$sigma))

  if (x$type == "active") {
    cat(sprintf("\nSequential testing results with alpha = %0.3f\n",x$alpha))
    tab = cbind(1:length(x$pv),x$vars,
      round(x$sign*x$vmat%*%x$y,3),
      round(x$sign*x$vmat%*%x$y/(x$sigma*sqrt(rowSums(x$vmat^2))),3),
      round(x$pv,3),round(x$ci,3))
    colnames(tab) = c("Step", "Var", "Coef", "Z-score", "P-value",
              "LowConfPt", "UpConfPt")
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


plot.fs <- function(x, breaks=TRUE, omit.zeros=TRUE, var.labels=TRUE, ...) {
  if (x$completepath) {
    k = length(x$action)+1
    beta = cbind(x$beta,x$bls)
  } else {
    k = length(x$action)
    beta = x$beta
  }
  p = nrow(beta)

  xx = 1:k
  xlab = "Step"

 if (omit.zeros) {
   good.inds = matrix(FALSE,p,k)
   good.inds[beta!=0] = TRUE
   changes = t(apply(beta,1,diff))!=0
   good.inds[cbind(changes,rep(F,p))] = TRUE
   good.inds[cbind(rep(F,p),changes)] = TRUE
   beta[!good.inds] = NA
  }

  plot(c(),c(),xlim=range(xx,na.rm=T),ylim=range(beta,na.rm=T),
       xlab=xlab,ylab="Coefficients",main="Forward stepwise path",...)
  abline(h=0,lwd=2)
  matplot(xx,t(beta),type="l",lty=1,add=TRUE)
  if (breaks) abline(v=xx,lty=2)
  if (var.labels) axis(4,at=beta[,k],labels=1:p,cex=0.8,adj=0)
  invisible()
}

