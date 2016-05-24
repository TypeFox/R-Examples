"rcnorm" <- function(n){
  sqrt(0.5)*(rnorm(n) + 1i*rnorm(n))
}

"dcmvnorm" <-
function (z, mean, sigma, log = FALSE) 
{
    if (is.vector(z)) {
        z <- matrix(z, ncol = length(z))
    }
    if (missing(mean)) {
        mean <- rep(0, length = ncol(z))
    }
    if (missing(sigma)) {
        sigma <- diag(ncol(z))
    }
    if (NCOL(z) != NCOL(sigma)) {
        stop("x and sigma have non-conforming size")
    }
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps))) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != NROW(sigma)) {
        stop("mean and sigma have non-conforming size")
    }

    distval <- quad.tdiag(solve(sigma),sweep(z,2,mean))
    
#    logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
    logdet <- sum(log(svd(z)$d))
    logretval <- Re(-(ncol(z) * log(pi) + logdet + distval))
    if (log) {
        return(logretval)
    } else {
        return(exp(logretval))
    }
}

"isHermitian" <-
function(x, tol= 100 * .Machine$double.eps){
  if(!isSymmetric(Re(x),tol=tol)){
    return(FALSE)
  }
  x <- Im(x)
  if(all(abs(x+t(x))<tol)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

"ishpd" <-
function(x,tol= 100 * .Machine$double.eps){ # "iphd" == Is Hermitian Positive Definite
  if(isHermitian(x) && all(zapim(svd(x)$d)>0)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

"rcmvnorm" <-
function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
          method = c("svd", "eigen", "chol"),
          tol= 100 * .Machine$double.eps){

  stopifnot(ishpd(sigma, tol)) # thus sigma known to be HPD
                                        
  if (length(mean) != nrow(sigma)) {
    stop("mean and sigma have non-conforming size")
  }
  
  method <- match.arg(method)
  if (method == "eigen") {
    ev <- eigen(sigma, symmetric=TRUE)
    D <- diag(sqrt(ev$values), ncol=length(ev$values))
    retval <- quad.tform(D, ev$vectors)
  } else if (method == "svd") {
    jj <- svd(sigma)
    if (!all(jj$d >= -sqrt(.Machine$double.eps) * abs(jj$d[1]))) {
      warning("sigma is numerically not positive definite")
    }
    D <- diag(sqrt(jj$d), ncol=length(jj$d))
    retval <- quad.tform(D, jj$u)
  } else {
    stop("option not recognized")
  }

  out <- matrix(rcnorm(n*ncol(sigma)),ncol=n)
  out <- Conj(crossprod(out,retval))
  out <- sweep(out, 2, mean, "+")
  colnames(out) <- names(mean)
  return(out)
}

"zapim" <-
function(x,tol= 100 * .Machine$double.eps){
  jj <- abs(Im(x))<tol
  if(all(jj)){
    return(Re(x))
  } else {
    Im(x[jj]) <- 0
    return(x)
  }
}

"Im<-" <-
function (x, value) 
{
    if (is.complex(value)) {
        stop("RHS must be pure real")
    }
    if (all(value == 0)) {
        return(Re(x))
    }
    else {
        return(Re(x) + (1i) * value)
    }
}

"Re<-" <-
function (x, value) 
{
    if (is.complex(value)) {
        stop("RHS must be pure real")
    }
    return((1i) * Im(x) + value)
}

"complex_CF" <- function(z1,z2,means,pos.def.matrix){
    exp(
        +0i
        +1i*Re(cprod(means,z2-z1))
        -quad.form(pos.def.matrix, z2-z1)
        )
}

"corr_complex" <-
function (z1, z2=NULL, distance.function = complex_CF, means=NULL, scales=NULL, pos.def.matrix=NULL){
    if(is.null(z2)){z2 <- z1}
    
    if (is.null(scales) & is.null(pos.def.matrix)) {
        stop("need either scales or a pos.definite.matrix")
    }
    if (!is.null(scales) & !is.null(pos.def.matrix)) {
        stop("scales *and* pos.def.matrix supplied.  corr() needs one only.")
    }
    if (is.null(pos.def.matrix)) {
        pos.def.matrix <- diag(scales, nrow = length(scales))
    }
    
    if(is.null(means)){
        means <- diag(pos.def.matrix)*0
    }
    
    out <-
        apply(z1, 1, function(y) {
            apply(z2, 1, function(x) {
                distance.function(x, y, means=means, pos.def.matrix=pos.def.matrix)
            })
        })
    
    return(as.matrix(out))
}

"interpolant.quick.complex" <-
function (x, d, zold, Ainv, scales = NULL, pos.def.matrix = NULL,  means=NULL,
    func = regressor.basis, give.Z = FALSE, distance.function = corr_complex, 
    ...) {
    if (is.null(scales) & is.null(pos.def.matrix)) {
        stop("need either scales or a pos.definite.matrix")
    }
    if (!is.null(scales) & !is.null(pos.def.matrix)) {
        stop("scales *and* pos.def.matrix supplied.  corr() needs one only.")
    }
    if (is.null(pos.def.matrix)) {
        pos.def.matrix <- diag(scales, nrow = length(scales))
    }
    x <- rbind(x)
    betahat <- betahat.fun(zold, Ainv, d, func = func)
    H <- regressor.multi(zold, func = func)
    mstar.star <- rep(NA, nrow(x))
    prior <- rep(NA, nrow(x))
    if (give.Z) {
        Z <- rep(NA, nrow(x))
        sigmahat.square <- sigmahatsquared(H, Ainv, d)
    }
    for (i in 1:nrow(x)) {
        hx <- func(x[i, ])
#        tx <- as.vector(corr_complex(x1 = zold, x2 = x[i, , drop = FALSE],
#                                     pos.def.matrix = pos.def.matrix,
#                                    means=means))
        tx <- as.vector(corr_complex(z1 = x[i, , drop = FALSE], z2=zold,
                                     pos.def.matrix = pos.def.matrix,
                                    means=means))
        prior[i] <- crossprod(hx, betahat)   # sic [that is, use crossprod() rather than cprod(), because hx is a row vector

        mstar.star[i] <- prior[i] + cprod(cprod(Ainv, tx), (d - H %*% betahat))
        if (give.Z) {
            cstar.x.x <- 1 - quad.form(Ainv, tx)
            cstar.star <- cstar.x.x + quad.form.inv(quad.form(Ainv, 
                H), hx - Conj(cprod(H, cprod(Ainv, tx))))
            Z[i] <- sqrt(abs(sigmahat.square * cstar.star))
        }
    }
    if (give.Z) {
        return(list(mstar.star = mstar.star, Z = Z, prior = prior))
    }
    else {
        return(mstar.star)
    }
}

# links to scales.likelihood() of the emulator package:
"scales.likelihood.complex" <-
function(pos.def.matrix, scales, means,  zold, z,  give_log = TRUE, func = regressor.basis){
    if(missing(pos.def.matrix)){
        pos.def.matrix <- diag(scales, nrow = length(scales))
    }
    H <- regressor.multi(zold, func = func)
    q <- ncol(H)
    n <- nrow(H)
    A <- corr_complex(z1=zold, means=means, pos.def.matrix=pos.def.matrix)
    Ainv <- solve(A)
    f <- function(M) {
        (-1) * sum(log(svd(M)$d))   # -0.5 for real case, -1 for complex
    }


#   bit1 <- log(sigmahatsquared(H, Ainv, d)) * (-(n - q)/2)   # real case
    bit1 <- log(sigmahatsquared(H, Ainv, z)) * (-(n - q)  )   # complex
    bit2 <- f(A)
    bit3 <- f(quad.form(Ainv, H))
    
    out <- drop(bit1 + bit2 + bit3)
    if (give_log) {
        return(out)
    }
    else {
        return(exp(out))
    }
}

"sd" <- function(x, na.rm = FALSE) {
  UseMethod("sd", x)  # thanks go to Kirill Mueller for advice
}

sd.complex <- function(x, na.rm = FALSE) {
    sqrt(drop(var(c(x),na.rm=na.rm)))
}


sd.default  <- stats::sd
var.default <- stats::var

"var" <- function(x,y=NULL,na.rm=FALSE,use){
    UseMethod("var",x)
}

"var.complex" <- function(x, y=NULL, na.rm = FALSE, use){
    if(na.rm){x <- x[!is.na(rowSums(x)),]}
    
    x <- as.matrix(x)
    x <- sweep(x,2,colMeans(x))
    if(is.null(y)){
       return(drop(zapim(cprod(x))/(nrow(x)-1)))
    } else {
       y <- as.matrix(y)
       if(na.rm){y <- y[!is.na(rowSums(y)),]}
       y <- sweep(y,2,colMeans(y))
       return(drop(zapim(cprod(x,y))/(nrow(x)-1)))
    }
}
