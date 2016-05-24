
## =============================================================================
## Solves a banded system of linear equations Ax=B
## =============================================================================

Solve.banded <- function(abd, nup, nlow,
             B=rep(0,times=ncol(abd)), full=(nrow(abd)==ncol(abd))) {

  B  <- as.matrix(B)
  Nb <- NCOL(B)

  Nmx  <- ncol(abd)
  nr   <- nrow(abd)

  if (full)   {    # full matrix was specified

    A <- abd
    if (nrow(abd) != ncol(abd))
      stop("cannot solve banded problem - nrows and ncols of abd are not the same, while the input matrix is said to be full")
    Aext <- rbind(matrix(data=0,nrow=nup,ncol=ncol(A)),
               A,matrix(data=0,nrow=nlow,ncol=ncol(A)))
    abd  <- matrix(nrow=nup+nlow+1,ncol=Nmx,data=
              Aext[(col(Aext))<=row(Aext)&col(Aext)>=row(Aext)-nlow-nup])
    nr   <- nrow(abd)
  }


  Nabd <- 2*nlow+nup+1
  if (nr != nlow+nup+1)
    stop("cannot solve banded problem - abd not compatible with nup and nlow")
  if (NROW(B) != Nmx)
    stop("cannot solve banded problem - B and abd not compatible")

  ABD <- rbind(matrix(data=0,nrow=nlow,ncol=Nmx),abd)
  IsError <- FALSE

    sol <-.Fortran("dgbsv",Nmx,as.integer(nlow),as.integer(nup),as.integer(Nb),
                   AB=ABD,LDAB=as.integer(nrow(ABD)),
                    ipiv=as.integer(rep(0,Nmx)),B=B,LDB=Nmx,info=as.integer(0))
  return(sol$B)

}
