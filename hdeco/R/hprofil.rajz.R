"hprofil.rajz" <-
function (BE=.HPROFIL, e=0.4, color1="lightgray", color2="darkgray") {

  NR <- dim(BE)[1]
  X <- 1:NR
  A <- BE[,1]
  BB <- BE[,2]
  B <- BE[,3]
  AMI <- min(A)
  AMA <- max(B+BB)
  
  # PREPARE AXES OF PLOT
  if(dim(.QKEP)[1]!=1){
    plot(x=c(-0.5, rep((X-1),2),NR-.5), y=c(0,A,BB+B,1), type="n", xlab="Step", ylab="Entropy", main="Entropy profile", ylim=c(AMI,AMA))
  }
  else{
    plot(x=c(-0.5,rep(X,2),NR+.5), y=c(0,A,BB+B,1), type="n", xlab="Step", ylab="Entropy", main="Entropy profile", ylim=c(AMI,AMA))
  }
  
  # DRAW PLOT BARS
  for (i in X){
    K <- BE[i,]
    if(dim(.QKEP)[1]!=1){
      polygon(c(i-1, i-e-1, i-e-1, i-1), c(0, 0, K[1], K[1]), col=color1)
      #polygon(c(i-1, i+e-1, i+e-1, i-1), c(K[2], K[2], K[2] + K[3], K[2] + K[3]), col=color2)
      polygon(c(i-1, i+e-1, i+e-1, i-1), c(0, 0, (K[2] + K[3]), (K[2] + K[3])), col=color2)
    }
    else {
      polygon(c(i,i-e,i-e,i), c(0,0,K[1],K[1]), col=color1)
      #polygon(c(i,i+e,i+e,i), c(K[2],K[2],K[2]+K[3],K[2]+K[3]), col=color2)
      polygon(c(i,i+e,i+e,i), c(0,0,K[2]+K[3],K[2]+K[3]), col=color2)
    }
  }
}

