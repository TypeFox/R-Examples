RWZAccept <-
function(a, first, last, constrained, impulses){
  QR <- qr(a)
  R <- qr.R(QR)
  Q <- qr.Q(QR)
  #
  if(R<0){Q  <- -1.0 * Q}
  #
  for(k in first:last){
    ik <- impulses[k, , ]%*%Q
    for(i in 1:length(constrained)){
      if(constrained[i]<0){
        value <- ik[-1.0*constrained[i]]
      }else{
        value <- -1.0 * ik[constrained[i]]
      }
      #
      if(value<0.0){
        if(k==first & i==1){
          Q <- 1.0 * Q
          ik <- 1.0 * ik
        } # comment this one out and uncomment bracket below.
      }else{
        acc <- 0
        rwz <- list(Q=Q, acc=acc, ika=ik)
        return(rwz)
        # } # comment out this
      }
    }
  }
  acc <- 1
  rwz <- list(Q=Q, acc=acc, ika=ik)
  return(rwz)
}
