"zsz" <-
function (ik=2,xremargi=F) {

  # WHAT DOES xremargi MEAN???

  MODEL <- sort(union(.AHIPO[[ik]],.NHIPO[[ik]]))
  QND <- .QND[MODEL]
  H <- length(QND)
  NAMES <- names(QND)

  if(xremargi) {
    XEK <- "X"==substring(NAMES,1,1)
    if(sum(XEK)==0)return(NULL)
    MODEL<-MODEL[!XEK]
    QND <- .QND[MODEL]
    H <- length(QND)
    NAMES <- names(QND)
  }

  QKEP <- bemar(MODEL)
  ZEK <- "Z"==substring(NAMES,1,1)

  if(sum(ZEK)==0) {
    return(NULL)
  }

  if(sum(!ZEK)==0) {
    return(NULL)
  }

  KI <- kifeltent(QKEP,QND,ZEK)
  if(is.null(KI)) return()
  QMARG <- bemar((1:H)[!ZEK],QKEP)
  ENTROPMARG <- entro(QMARG)
  NR <- dim(KI)[1]
  KI <- cbind(rep(ENTROPMARG,NR),KI)
  dimnames(KI)[[2]][1] <- "HMARG"
  KI
  cat("\nLeaving zsz()")
}

