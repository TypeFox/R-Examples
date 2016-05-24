"rot.lsq" <-function(xx,
                    yy,
                    xfit=rep(TRUE,length(xx)),
                    yfit=xfit,
                    verbose=FALSE) {
  
  # Coordinate superposition with the Kabsch algorithm
  # from Acta Cryst (1978) A34 pp827-828 (to which equation no. refer)
  # yy is the target (i.e. fixed)

  xx <- matrix(xx,nrow=3, ) 
  x <- matrix(xx[xfit],nrow=3, ) 
  y <- matrix(yy[yfit],nrow=3, )
  
  if(length(x) != length(y)) stop("dimension mismatch in x and y")

  # mean positions 
  xbar <- apply(x,1,mean) ; ybar <- apply(y,1,mean)

  # center both sets
  xx <- sweep(xx,1,xbar) # NB xx centred on xbar 
  x <- sweep(x,1,xbar) ; y <- sweep(y,1,ybar)
  
  #irmsd <- sqrt(sum((x-y)^2)/dim(y)[2])
  #cat("#irmsd= ",round(irmsd,6),"\n")

  # generate the 3x3 moment matrix: R (Equation 3)
  R <- y %*% t(x)

  # form R'R
  RR <- t(R) %*% R
  
  # diagonalize R'R
  prj <- eigen(RR)
  prj$values[prj$values < 0 & prj$values >= -1.0E-12]<-1.0E-12
  
  # form A
  A <- prj$vectors

  # make explicitly rh system
  #	A[,3] <- v3cross(A[,1],A[,2])
  # inline the cross-product function call.
  b<-A[,1]; c <- A[,2]
  A[1,3] <- (b[2] * c[3]) - (b[3] * c[2])
  A[2,3] <- (b[3] * c[1]) - (b[1] * c[3])
  A[3,3] <- (b[1] * c[2]) - (b[2] * c[1])

  # form B (==RA) (Equation 8)
  B <- R %*% A

  # normalize B
  # B <- sweep(B,2,sqrt(apply(B^2,2,sum)),"/")
  B <- sweep(B,2,sqrt(prj$values),"/")

  # make explicitly rh system
  #	B[,3] <- v3cross(B[,1],B[,2])
  # inline the cross-product function call.
  b<-B[,1]; c <- B[,2]
  B[1,3] <- (b[2] * c[3]) - (b[3] * c[2])
  B[2,3] <- (b[3] * c[1]) - (b[1] * c[3])
  B[3,3] <- (b[1] * c[2]) - (b[2] * c[1])

  # form U (==Ba) (Equation 7)
  # U is the rotation matrix
  U <-  B %*% t(A)

  # here we apply transformation matrix to *all* elements of xx
  # rotate xx (Uxx)
  xx <- U %*% xx

  if(verbose) {
    ## also apply it to the subset, in order to compute residual
    x <- U %*% x

    ## estimate of residuals 
    frmsd <- sqrt(sum((x-y)^2)/dim(y)[2])
    cat("#rmsd= ",round(frmsd,6),"\n")
  }
# fest <- iest - sum(sqrt(prj$values))
# print(sqrt((2*fest)/dim(y)[2]))

# return xx centred on y
  xx <- sweep(xx,1,ybar,"+")
  as.vector(xx)
}

