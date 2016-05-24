BPEC.circles <- function(x=0, y=0, radius=1, col='blue', border=par("bg"),angle=360)
  {
    theta <- 2* pi * seq(0, angle, by=20) / 360
    ctheta <- cos(theta)
    stheta <- sin(theta)
    m<-c(x,x+ctheta*radius)
    n<-c(y,y+stheta*radius)
    polygon(m,n,col=col, border=border)
}
