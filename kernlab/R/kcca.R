## Simple kernel canonical corelation analysis
## author: alexandros karatzoglou

setGeneric("kcca",function(x, y, kernel="rbfdot",  kpar=list(sigma = 0.1), gamma=0.1, ncomps = 10, ...) standardGeneric("kcca"))
setMethod("kcca", signature(x = "matrix"),
          function(x,y,kernel="rbfdot",kpar=list(sigma=0.1), gamma=0.1, ncomps =10, ...)
          {
            x <- as.matrix(x)
            y <- as.matrix(y)

            if(!(nrow(x)==nrow(y)))
                stop("Number of rows in x, y matrixes is not equal")
            if(!is(kernel,"kernel"))
              {
                if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
                kernel <- do.call(kernel, kpar)
              }
            if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

            Kx <- kernelMatrix(kernel,x)
            Ky <- kernelMatrix(kernel,y)
            
            n <- dim(Kx)[1]
            m <- 2
            ## Generate LH
            VK <- matrix(0,n*2,n);

            VK[0:n,] <- Kx
            VK[(n+1):(2*n),] <- Ky
            LH <- tcrossprod(VK, VK)

            for (i in 1:m)
              LH[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)] <- 0

            ## Generate RH
            RH <- matrix(0,n*m,n*m)
            RH[1:n,1:n] <- (Kx + diag(rep(gamma,n)))%*%Kx + diag(rep(1e-6,n))
            RH[(n+1):(2*n),(n+1):(2*n)] <- (Ky + diag(rep(gamma,n)))%*%Ky + diag(rep(1e-6,n))
            RH <- (RH+t(RH))/2

            ei <- .gevd(LH,RH)

            ret <- new("kcca")

            kcor(ret) <- as.double(ei$gvalues[1:ncomps])
            xcoef(ret) <- matrix(as.double(ei$gvectors[1:n,1:ncomps]),n)
            ycoef(ret) <- matrix(as.double(ei$gvectors[(n+1):(2*n),1:ncomps]),n)
            ## xvar(ret) <- rotated(xpca) %*% cca$xcoef
            ## yvar(ret) <- rotated(ypca) %*% cca$ycoef
            return(ret)
          })

## gevd compute the generalized eigenvalue 
## decomposition for (a,b)
.gevd<-function(a,b=diag(nrow(a))) {
  bs<-.mfunc(b,function(x) .ginvx(sqrt(x)))
  ev<-eigen(bs%*%a%*%bs)
  return(list(gvalues=ev$values,gvectors=bs%*%ev$vectors))
}

## mfunc is a helper to compute matrix functions
.mfunc<-function(a,fn=sqrt) {
  e<-eigen(a); y<-e$vectors; v<-e$values
  return(tcrossprod(y%*%diag(fn(v)),y))
}

## ginvx is a helper to compute reciprocals
.ginvx<-function(x) {ifelse(x==0,0,1/x)}

