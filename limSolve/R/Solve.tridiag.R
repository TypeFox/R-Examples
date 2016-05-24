
## =============================================================================
## Solves tridiagonal system of linear equations
## =============================================================================

Solve.tridiag  <- function(diam1, dia, diap1,
                  B=rep(0,times=length(dia)))  {

  B  <- as.matrix(B)
  Nb <- NCOL(B)

  Nmx <- length(dia)
  if (length(diam1) != Nmx-1)
    stop("cannot solve tridiagonal problem - diam1 and dia not compatible")
  if (length(diap1) != Nmx-1)
    stop("cannot solve tridiagonal problem - diap1 and dia not compatible")
  if (NROW(B  ) != Nmx)
    stop("cannot solve tridiagonal problem - B and dia not compatible")
  value=rep(0.,Nmx)
#  sol <-.Fortran("tridia",Nmx=Nmx,
#         au=as.double(c(0,diam1)),bu=as.double(dia),cu=as.double(c(diap1,0)),
#         du=as.double(B),value=as.double(value))

#  return(sol$value)    # vector with solution of tridiagonal system


#DGTSV( N, NB, DL, D, DU, B, LDB, INFO )
  sol <-.Fortran("dgtsv",N=Nmx,nrhs=Nb,DL=as.double(diam1),
        D=as.double(dia),DU=as.double(diap1),B=as.double(B),LDB=Nmx,
        INFO=as.integer(0) )
  result <- matrix(nrow = nrow(B), data = sol$B)
  return(result)
}
