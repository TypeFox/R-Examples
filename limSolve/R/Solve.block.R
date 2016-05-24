Solve.block <- function(Top,AR,Bot,B,overlap) {

  B  <- as.matrix(B)
  cB <- NCOL(B)
  rB <- NROW(B)
  BB <- B[,1]

  NrowTop <- nrow(Top)
  NrowBot <- nrow(Bot)
  Dim <- dim(AR)
  NrwBlk <- Dim[1]
  NclBlk <- Dim[2]
  Nbloks <- Dim[3]
  N <- Nbloks*NrwBlk + overlap
  if (N != rB)
    stop(paste("AR and B not compatible: nrow B should be: ",N))
  X <- matrix(data=as.double(0.), nrow=N, ncol=cB)


  sol <- .Fortran("block", N=as.integer(N), TOP=Top, NRWTOP=NrowTop,
    NOVRLP=as.integer(overlap),AR=AR,NRWBLK=NrwBlk,NCLBLK=NclBlk,NBLOKS=Nbloks,
    BOT=Bot,NRWBOT=NrowBot,PIVOT=as.integer(rep(0,N)),rB=rB,
    cB=cB, B=as.double(B), X=X, IFLAG=as.integer(0),
    Tmp1=rep(0.,rB),Tmp2=rep(0.,N))

  return(sol$X)
}



