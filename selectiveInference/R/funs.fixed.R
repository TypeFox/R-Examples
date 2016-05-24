# Lasso inference function (for fixed lambda). Note: here we are providing inference
# for the solution of
# min 1/2 || y - \beta_0 - X \beta ||_2^2 + \lambda || \beta ||_1

fixedLassoInf <- function(x, y, beta, lambda, intercept=TRUE, sigma=NULL, alpha=0.1,
                     type=c("partial","full"), tol.beta=1e-5, tol.kkt=0.1,
                     gridrange=c(-100,100), bits=NULL, verbose=FALSE) {
  
  this.call = match.call()
  type = match.arg(type)
  checkargs.xy(x,y)
  if (missing(beta) || is.null(beta)) stop("Must supply the solution beta")
  if (missing(lambda) || is.null(lambda)) stop("Must supply the tuning parameter value lambda") 
  checkargs.misc(beta=beta,lambda=lambda,sigma=sigma,alpha=alpha,
                 gridrange=gridrange,tol.beta=tol.beta,tol.kkt=tol.kkt)
  if (!is.null(bits) && !requireNamespace("Rmpfr",quietly=TRUE)) {
    warning("Package Rmpfr is not installed, reverting to standard precision")
    bits = NULL
  }
  
  n = nrow(x)
  p = ncol(x)
  beta = as.numeric(beta)
  if (length(beta) != p) stop("beta must have length equal to ncol(x)")

  # If glmnet was run with an intercept term, center x and y
  if (intercept==TRUE) {
    obj = standardize(x,y,TRUE,FALSE)
    x = obj$x
    y = obj$y
  }

  # Check the KKT conditions
  g = t(x)%*%(y-x%*%beta) / lambda
  if (any(abs(g) > 1+tol.kkt * sqrt(sum(y^2))))
    warning(paste("Solution beta does not satisfy the KKT conditions",
                  "(to within specified tolerances)"))

  vars = which(abs(beta) > tol.beta / sqrt(colSums(x^2)))
  if(length(vars)==0){
      cat("Empty model",fill=T)
      return()
  }
  if (any(sign(g[vars]) != sign(beta[vars])))
    warning(paste("Solution beta does not satisfy the KKT conditions",
                  "(to within specified tolerances). You might try rerunning",
                  "glmnet with a lower setting of the",
                  "'thresh' parameter, for a more accurate convergence."))
  
  # Get lasso polyhedral region, of form Gy >= u
  out = fixedLasso.poly(x,y,beta,lambda,vars)
  G = out$G
  u = out$u

  # Check polyhedral region
  tol.poly = 0.01 
  if (min(G %*% y - u) < -tol.poly * sqrt(sum(y^2)))
    stop(paste("Polyhedral constraints not satisfied; you must recompute beta",
               "more accurately. With glmnet, make sure to use exact=TRUE in coef(),",
               "and check whether the specified value of lambda is too small",
               "(beyond the grid of values visited by glmnet).",
               "You might also try rerunning glmnet with a lower setting of the",
               "'thresh' parameter, for a more accurate convergence."))

  # Estimate sigma
  if (is.null(sigma)) {
    if (n >= 2*p) {
      oo = intercept
      sigma = sqrt(sum(lsfit(x,y,intercept=oo)$res^2)/(n-p-oo))
    }
    else {
      sigma = sd(y)
      warning(paste(sprintf("p > n/2, and sd(y) = %0.3f used as an estimate of sigma;",sigma),
                    "you may want to use the estimateSigma function"))
    }
  }
 
  k = length(vars)
  pv = vlo = vup = numeric(k) 
  vmat = matrix(0,k,n)
  ci = tailarea = matrix(0,k,2)
  sign = numeric(k)

  if (type=="full" & p > n)
      warning(paste("type='full' does not make sense when p > n;",
                    "switching to type='partial'"))
  
  if (type=="partial" || p > n) {
    xa = x[,vars,drop=F]
    M = pinv(crossprod(xa)) %*% t(xa)
  }
  else {
    M = pinv(crossprod(x)) %*% t(x)
    M = M[vars,,drop=F]
  }
  
  for (j in 1:k) {
    if (verbose) cat(sprintf("Inference for variable %i ...\n",vars[j]))
    
    vj = M[j,]
    mj = sqrt(sum(vj^2))
    vj = vj / mj        # Standardize (divide by norm of vj)
    sign[j] = sign(sum(vj*y))
    vj = sign[j] * vj
    a = poly.pval(y,G,u,vj,sigma,bits)
    pv[j] = a$pv 
    vlo[j] = a$vlo * mj # Unstandardize (mult by norm of vj)
    vup[j] = a$vup * mj # Unstandardize (mult by norm of vj)
    vmat[j,] = vj * mj  # Unstandardize (mult by norm of vj)

    a = poly.int(y,G,u,vj,sigma,alpha,gridrange=gridrange,
      flip=(sign[j]==-1),bits=bits)
    ci[j,] = a$int * mj # Unstandardize (mult by norm of vj)
    tailarea[j,] = a$tailarea
  }
  
  out = list(type=type,lambda=lambda,pv=pv,ci=ci,
    tailarea=tailarea,vlo=vlo,vup=vup,vmat=vmat,y=y,
    vars=vars,sign=sign,sigma=sigma,alpha=alpha,
    call=this.call)
  class(out) = "fixedLassoInf"
  return(out)
}

#############################


fixedLasso.poly=
function(x, y, beta, lambda, a) {
  xa = x[,a,drop=F]
  xac = x[,!a,drop=F]
  xai = pinv(crossprod(xa))
  xap = xai %*% t(xa)
  za = sign(beta[a])
  if (length(za)>1) dz = diag(za)
  if (length(za)==1) dz = matrix(za,1,1)

  P = diag(1,nrow(xa)) - xa %*% xap
  #NOTE: inactive constraints not needed below!
 
  G = -rbind(
   #   1/lambda * t(xac) %*% P,
   # -1/lambda * t(xac) %*% P,
    -dz %*% xap
      )
     lambda2=lambda
     if(length(lambda)>1) lambda2=lambda[a]
  u = -c(
   #   1 - t(xac) %*% t(xap) %*% za,
   #   1 + t(xac) %*% t(xap) %*% za,
    -lambda2 * dz %*% xai %*% za)

  return(list(G=G,u=u))
}
# Moore-Penrose pseudo inverse for symmetric matrices

pinv <- function(A, tol=.Machine$double.eps) {
  e = eigen(A)
  v = Re(e$vec)
  d = Re(e$val)
  d[d > tol] = 1/d[d > tol]
  d[d < tol] = 0
  if (length(d)==1) return(v*d*v)
  else return(v %*% diag(d) %*% t(v))
}

##############################

print.fixedLassoInf <- function(x, tailarea=TRUE, ...) {
  cat("\nCall:\n")
  dput(x$call)

  cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n",
              x$sigma))
  
  cat(sprintf("\nTesting results at lambda = %0.3f, with alpha = %0.3f\n",x$lambda,x$alpha))
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
 
  cat(sprintf("\nNote: coefficients shown are %s regression coefficients\n",
              ifelse(x$type=="partial","partial","full")))
  invisible()
}

#estimateLambda <- function(x, sigma, nsamp=1000){
#  checkargs.xy(x,rep(0,nrow(x)))
#  if(nsamp < 10) stop("More Monte Carlo samples required for estimation")
#  if (length(sigma)!=1) stop("sigma should be a number > 0")
 # if (sigma<=0) stop("sigma should be a number > 0")
                                      
 # n = nrow(x)
 # eps = sigma*matrix(rnorm(nsamp*n),n,nsamp)
 # lambda = 2*mean(apply(t(x)%*%eps,2,max))
 # return(lambda)
#}
    
