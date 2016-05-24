# ============================================ #
# Functions for solving least squares problems #
#                                              #
# Author: Yong Wang (yongwang@auckland.ac.nz)  #
#         Department of Statistics             #
#         University of Auckland               #
#         New Zealand                          #
# ============================================ #

# ----------------------------------------------------------------------
# This package provides the R interface functions to some Fortran
# subroutines in the LSEI package and implements some algorithms described
# in Lawson and Hanson (1974), "Solving Least Squares Problems",
# Prentice-HalL.  The source code of this package is downloaded from
# http://www.netlib.org/lawson-hanson/all

# Last modified: Mon Sep 21 11:43:14 NZST 2015
# ----------------------------------------------------------------------

# .First.lib <- function(lib, pkg) {
#   library.dynam("lsei", pkg, lib)
# }

#---------------------------------- #
# Nonnegative least squares (NNLS): #
#                                   #
#        Minimize     ||ax - b||    #
#        subject to   x >= 0        # 
# --------------------------------- #

nnls = function(a, b) {
  m = as.integer(dim(a)[1])
  n = as.integer(dim(a)[2])
  storage.mode(a) = "double"
  storage.mode(b) = "double"
  x = double(n)                       # only for output
  rnorm = double(1)                   # only for output
  w = x                               # n-vector of working space
  zz = b                              # m-vector of working space
  index = integer(n)                  # n-vector index, only for output
  mode = integer(1)                   # success-failure flag; = 1, success
  .Fortran("nnls",a,m,m,n,b,x=x,rnorm=rnorm,w,zz,index=index,
           mode=mode,PACKAGE="lsei")[c("x","rnorm","index","mode")]
}

# Solving a partial NNLS problem

# Input:
#
# a      Matrix
# b      Vector or a one-column matrix
# k      The first k variables are not restricted to non-negativity. 

# Output:
#
# x      Solution vector (in the original order of variables)
# r      Upper triangular matrix
# b      Vector b after Q transformation
# rnorm  Euclidean norm of the residual vector
# mode   Termination mode:
#           = 1, successful;
#           = 2, bad dimensions of the problem, either m <= 0 or n <= 0;
#           = 3, iteration count exceeded (more than 3*n iterations).
# k      Number of free variables after eliminating linear dependence
# sum    = NULL, if NN-restricted coefficients are not further
#            restricted to have a fixed sum;
#        = positive value, if NN-restricted coefficients are further
#            restricted to have a fixed positive sum;

pnnls = function(a, b, k=0, sum=NULL) {
  if(!is.vector(b)) b = drop(b)
  m = as.integer(dim(a)[1])
  n = as.integer(dim(a)[2])
  if(!is.null(sum)) {
    if(sum <= 0) stop("Argument 'sum' must be positive or NULL")
    if(k<n) a[,(k+1):n] = a[,(k+1):n] * sum - b
    else stop("k == ncol(a) (simplex contains no variable)")
    a = rbind(a, c(rep(0,k), rep(1, n-k)))
    b = c(rep(0, m), 1)
    m = as.integer(m+1)
  }
  storage.mode(a) = "double"
  storage.mode(b) = "double"
  x = double(n)                       # only for output
  rnorm = double(1)                   # only for output
  w = x                               # n-vector of working space
  zz = b                              # m-vector of working space
  index = integer(n)                  # n-vector index, only for output
  mode = integer(1)                   # success-failure flag; = 1, success
  k = as.integer(k)
  r = .Fortran("pnnls",r=a,m,m,n,b=b,x=x,rnorm=rnorm,w,zz,index=index,
      mode=mode,k=k,PACKAGE="lsei")
  r$r = r$r[1:min(m,n),]
  if(!is.null(sum)) {
    r$x = r$x / sum(r$x[(r$k+1):n])
    r$x[(r$k+1):n] = r$x[(r$k+1):n] * sum
  }
  r[c("x","r","b","rnorm","index","mode","k")]
}

# --------------------------------- #
# Least distance programming (LDP): #
#                                   #
#        Minimize  ||x||            #
#        Suject to e x >= f         #
# --------------------------------- #

ldp = function(e, f, tol=1e-15) {
  k = ncol(e)           # number of variables
  E = rbind(t(e), t(f))
  h = c(rep(0, k), 1)
  r = E %*% nnls(E, h)$x - h
  if( sum(r^2) <= tol ) stop("Incompatible inequalities in ldp()")
  as.vector( -r[1:k] / r[k+1] )
}

# # Example from Lawson and Hanson. (1974), p.171:
# e = cbind(c(-.207,-.392,.599), c(2.558, -1.351, -1.206))
# f = c(-1.3,-.084,.384)
# ldp(e, f)           # Solution: 0.1268538 -0.2554018

# G = matrix(rnorm(12), nrow=4); h = rnorm(4); x = ldp(G, h); print(x); G %*% x - h

# --------------------------------------------------------------- #
# Least squares problem with linear inequality constraints (LSI): #
#                                                                 #
#         Minimize       || a x - b ||                            #
#         Subject to     e x >= f                                 #
# --------------------------------------------------------------- #

lsi = function(a, b, e=NULL, f=NULL, lower=-Inf, upper=Inf) {
  if(is.vector(e)) dim(e) = c(1, length(e))
  if(any(lower != -Inf, upper != Inf)) {
    k0 = ncol(a)
    lower = rep(lower, len=k0)
    upper = rep(upper, len=k0)
    jl = lower != -Inf
    ju = upper != Inf
    e = rbind(e, diag(1, k0)[jl,], diag(-1, k0)[ju,])
    f = c(f, lower[jl], -upper[ju])
  }
  if(is.null(e)) return( pnnls(a, b, k=ncol(a))$x )
  a.svd = svdrs(a, b)
  k = sum( a.svd$d > max(a.svd$d) * 1e-10 )     # pseudo-rank
  Q1b = a.svd$uTb[1:k,,drop=FALSE]
  K1 = a.svd$v[,1:k,drop=FALSE]
  et = e %*% K1 %*% diag( 1/a.svd$d[1:k], nrow=k )
  ft = f - et %*% Q1b
  drop( K1 %*% ((ldp(et, ft) + Q1b) / a.svd$d[1:k]) )
}

# # Example from Lawson and Hanson. (1974), p.170:
# a = cbind(c(.25,.5,.5,.8),rep(1,4)); b = c(.5,.6,.7,1.2); e = cbind(c(1,0,-1),c(0,1,-1)); f = c(0,0,-1); lsi(a, b, e, f)    # Solution: 0.6213152 0.3786848

# E = matrix(rnorm(24), nrow=8); f = rnorm(8); G = matrix(rnorm(12), nrow=4); h = rnorm(4); x = lsi(E, f, G, h); G %*% x - h

# -------------------------------------------------------------------------- #
# Least squares problem with both linear equalities and inequalities (LSEI): #
#                                                                            #
#          Minimize      || a x - b ||                                       #
#          Subject to    c x = d                                             #
#          and           e x >= f                                            #
# -------------------------------------------------------------------------- #

lsei = function(a, b, c=NULL, d=NULL, e=NULL, f=NULL, lower=-Inf, upper=Inf) {
  if(is.null(c)) return(lsi(a, b, e, f))
  if(is.vector(c)) dim(c) = c(1, length(c))
  k = nrow(c)
  c.qr = qr( t(c) )
  L = t(qr.R(c.qr))          # To get back c:   L %*% t(qr.Q(c.qr))
  y1 = forwardsolve(L, d)
  ad = t(qr.qty(c.qr, t(a)))
  if(any(lower != -Inf, upper != Inf)) {
    k0 = ncol(a)
    lower = rep(lower, len=k0)
    upper = rep(upper, len=k0)
    jl = lower != -Inf
    ju = upper != Inf
    e = rbind(e, diag(1, k0)[jl,], diag(-1, k0)[ju,])
    f = c(f, lower[jl], -upper[ju])
  }
  if(is.null(e)) 
    y2 = pnnls(ad[,-(1:k),drop=FALSE], b - ad[,1:k,drop=FALSE] %*% y1,
        ncol(ad)-k)$x
  else {
    if(is.vector(e)) dim(e) = c(1, length(e))
    ed = t(qr.qty(c.qr, t(e)))
    y2 = lsi(ad[,-(1:k),drop=FALSE], b - ad[,1:k,drop=FALSE] %*% y1,
        ed[,-(1:k),drop=FALSE], f - ed[,1:k,drop=FALSE] %*% y1)
  }
  qr.qy(c.qr, c(y1, y2))
}

# beta = c(rnorm(2), 1); beta[beta<0] = 0; beta = beta/sum(beta)
# a = matrix(rnorm(18), ncol=3); b = a %*% beta + rnorm(3,sd=.1); c = matrix(rep(1, 3), nrow=1); d = 1; e = diag(rep(1,3)); f = rep(0,3); lsei(a, b, c, d, e, f)

# # c = matrix(rnorm(6), ncol=3); d = rnorm(2); a = matrix(rnorm(24), nrow=8); b = rnorm(8); e = matrix(rnorm(12), nrow=4); f = rnorm(4); x = lsei(a, b, c, d, e, f); print(x); print(c %*% x - d); e %*% x - f

# ------------------------------------------------------- #
# Least squares solution using Householder transformation #
# ------------------------------------------------------- #

hfti = function(a, b, tol = 1e-7) {
  if ( is.vector(b) ) b = as.matrix(b)
  if ( !(is.matrix(a) & is.matrix(b)) ) stop("a or b not a matrix")
  m = as.integer(dim(a)[1])
  n = as.integer(dim(a)[2])
  if ( m !=  dim(b)[1] ) stop("dim(a)[1] != dim(b)[1]")
  nb = as.integer(dim(b)[2])
  storage.mode(a) = "double"
  storage.mode(b) = "double"
  krank = as.integer(0)
  rnorm = as.double(rep(0,nb))        # only for output
  h = g = as.double(rep(0,n))         # n-vector of working space
  ip = rep.int(0,n)                   # m-vector of working space
  .Fortran("hfti",a=a,m,m,n,b=b,m,nb,tol,krank=krank,rnorm=rnorm,h,g,
           ip=ip,PACKAGE="lsei")[c("b","krank","rnorm")]
}

# ---------------------------------------------------------- #
# svdrs: singular value decomposition with right side vector #
#                                                            #
# For the least squares problem                              #
#                                                            #
#                ||a x - b||                                 #
# ---------------------------------------------------------- #

svdrs = function(a, b) {
  m1 = nrow(a)
  n1 = ncol(a)
  if( m1 < n1 ) a = rbind(a, matrix(0, nrow=n1-m1, ncol=n1))
  mda = nrow(a)
  k = min(m1, n1)
  missing.b = FALSE
  if( missing(b) ) {b = diag( rep(1.0, mda), nrow=mda ); missing.b = TRUE}
  if( is.vector(b) ) b = as.matrix(b)
  nb = ncol(b)
  if( nrow(b) < mda ) b = rbind(b, matrix(0, nrow=mda-nrow(b), ncol=ncol(b)))
  s = double(n1)
  work = double(2*n1)
  storage.mode(a) = "double"
  storage.mode(b) = "double"
  r = .Fortran("svdrs",a=a,mda,m1,n1,b=b,nrow(b),nb,s=s,
      work,PACKAGE="lsei")[c("s","a","b")]
  if( missing.b )
    list(d=r$s[1:k], u=t(r$b[1:min(k,nrow(b)), 1:min(m1,nb),drop=FALSE]),
         v=r$a[1:n1,1:k,drop=FALSE])
  else list(d=r$s[1:k], uTb=r$b[1:min(k,nrow(b)), 1:min(m1,nb),drop=FALSE],
            v=r$a[1:n1,1:k,drop=FALSE])
}

# x = matrix(rnorm(6), nrow=3)
# r = svdrs(x)
# r$u %*% diag(r$d) %*% t(r$v) - x
# svdrs(x, 1:3)

# Quadratic programming

# p^T x + x^T q x / 2

# h is positive semidefinite

qp = function(q, p, c=NULL, d=NULL, e=NULL, f=NULL,
    lower=-Inf, upper=Inf, tol=1e-15) {
  eq = eigen(q)
  v2 = sqrt(eq$values[eq$values >= eq$values[1] * tol])
  kr = length(v2)
  a = t(eq$vectors[,1:kr,drop=FALSE]) * v2
  b = - colSums(eq$vectors[,1:kr,drop=FALSE] * p / rep(v2, each=length(p)))
  lsei(a, b, c, d, e, f, lower, upper)
}

# partial nonnegativity quadratic programming

pnnqp = function(q, p, k=0, sum=NULL, tol=1e-15) {
  eq = eigen(q)
  v2 = sqrt(eq$values[eq$values >= eq$values[1] * tol])
  kr = length(v2)
  a = t(eq$vectors[,1:kr,drop=FALSE]) * v2
  b = - colSums(eq$vectors[,1:kr,drop=FALSE] * p / rep(v2, each=length(p)))
  pnnls(a, b, k, sum)
}
