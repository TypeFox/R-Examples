buildScoreMatrix = function() {
  # Get all of the constants
  a = getConstants()
  R1 = a$R1 
  R2 = a$R2
  R3 = a$R3
  R4 = a$R4
  R5 = a$R5
  R = a$R 
  S = a$S
  
  out = .C("BuildScoreMatrix",
    A=as.integer(numeric(681*681)),
    R=as.double(c(R1,R2,R3,R4,R5,R)),
    S=as.integer(S),
    package="darts")

  return(matrix(out$A,nrow=681))  
}

simpleExpScores = function(s) {
  # Build the score matrix
  A = buildScoreMatrix()
  
  # Next build the Gaussian density (1d) matrices
  b1 = 1/sqrt(2*pi*s)*exp(-(-340:340)^2/(2*s))
  b2 = 1/sqrt(2*pi*s)*exp(-(-340:340)^2/(2*s))

  # Build C
  C = t(mvfft(mvfft(t(A))*fft(b2), inverse=TRUE)/681)
  C = Re(C[,c(341:681,1:340)])
 
  # Finally build the convolution
  E = mvfft(mvfft(C)*fft(b1), inverse=TRUE)/681
  E = Re(E[c(341:681,1:340),])

  return(E[171:511,171:511])
}

generalExpScores = function(Sig) {
  # Build the score matrix
  A = buildScoreMatrix()

  # Next build the Gaussian density matrix
  x = matrix(rep(-340:340,681),nrow=681)
  y = matrix(rep(-340:340,681),nrow=681,byrow=TRUE)
  det = Sig[1]*Sig[2]-Sig[3]^2
  B = 1/(2*pi*sqrt(det))*
    exp(-(Sig[2]*x^2 - 2*Sig[3]*x*y + Sig[1]*y^2)/(2*det))
  
  E = fft(fft(A)*fft(B),inverse=TRUE)/681^2
  E = Re(E[c(341:681,1:340),c(341:681,1:340)])

  return(E[171:511,171:511])
}

drawHeatmap = function(e,col=heat.colors(30)) {
  R = getConstants()$R
  par(mar=c(0,0,0,0))
  image(-R:R,-R:R,e,axes=FALSE,xlim=c(-R-30,R+30),ylim=c(-R-30,R+30),col=col)
  rect(-R,-R,R,R,lwd=1,border="black")
}

drawAimSpot = function(e,col="blue",pch=19,...) {
  R = getConstants()$R
  a = argmax(e)
  points(a[1]-R,a[2]-R,col=col,pch=pch,...)
}

argmax = function(e) {
  k = which.max(e)
  i = (k-1)%%nrow(e) + 1
  j = floor(k/nrow(e)) + 1 
  return(c(i,j))
}
