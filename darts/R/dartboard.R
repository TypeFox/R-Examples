# Draws the dartboard, either a new plot or on top of an existing one
drawBoard = function(new=FALSE, lines=TRUE, numbers=TRUE, outside=TRUE, col="black", ...) {
  # Get all of the constants
  a = getConstants()
  R1 = a$R1 
  R2 = a$R2
  R3 = a$R3
  R4 = a$R4
  R5 = a$R5
  R = a$R 
  S = a$S
  
  if (new) {
    par(mar=c(0,0,0,0))
    plot(c(),c(),axes=FALSE,xlim=c(-R-15,R+15),ylim=c(-R-15,R+15))
  }

  if (lines) {
    t = seq(0,2*pi,length=5000)
    x = cos(t)
    y = sin(t)
    points(R*x,R*y,col=col,type="l")
    points(R5*x,R5*y,col=col,type="l")
    points(R4*x,R4*y,col=col,type="l")
    points(R3*x,R3*y,col=col,type="l")
    points(R2*x,R2*y,col=col,type="l")
    points(R1*x,R1*y,col=col,type="l")
    
    t0 = pi/2 + 2*pi/40
    points(c(R2*cos(t0),R*cos(t0)), 
           c(R2*sin(t0),R*sin(t0)),col=col,type="l")
    for (i in 1:19) {
      t1 = t0 - i*2*pi/20
      points(c(R2*cos(t1),R*cos(t1)), 
             c(R2*sin(t1),R*sin(t1)),col=col,type="l")
    }
  }

  if (numbers) {
    if (outside) r = R+10
    else r = 0.8*R
    for (i in 1:20) {
      t1 = pi/2 - (i-1)*2*pi/20
      text(r*cos(t1),r*sin(t1),S[i],...)
    }
  }
}
 
