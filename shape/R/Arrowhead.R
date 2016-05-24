
##==============================================================================
## Arrowhead    : draws arrowhead, various shapes
##==============================================================================

Arrowhead <- function(x0, y0, angle=0, arr.length=0.4,
  arr.width=arr.length/2, arr.adj=0.5, arr.type="curved",
  lcol="black", lty=1, arr.col=lcol, arr.lwd = 2, npoint = 5) {

 ## points of polygon, as drawn in graph with x- and y- ranges -5,5

  if ( arr.type=="curved") { # composed as section of circels

    rad <- 0.7                                        # radius of outer circles
    len <- 0.25*pi
    mid <- c(0,rad)

    x   <- seq(1.5*pi+len,1.5*pi,length.out=npoint)
    rr  <- cbind(mid[1]-rad*cos(x),mid[2]+rad*sin(x)) #part of circle
    mid <- c(0,-rad)
    x   <- rev(x)
    rr  <- rbind(rr,cbind(mid[1]-rad*cos(x),mid[2]-rad*sin(x)))
    mid <-c(rr[nrow(rr),1],0)
    rd  <-rr[1,2]
    x   <-seq(pi/2,3*pi/2,length.out=3*npoint)        #part of ellipse
    rr <- rbind(rr,cbind(mid[1]-rd*0.25*cos(x),mid[2]-rd*sin(x)))
    rr[,1] <- rr[,1]*2.6
    rr[,2] <- rr[,2]*3.45
  } else

  if (arr.type=="triangle") {
    x   <- c(-0.2,0.0,-0.2)
    y   <- c(-0.1,0.0,0.1)
    rr  <- 6.22*cbind(x,y)
  } else

  if (arr.type %in% c("circle","ellipse") )  {

    if (arr.type=="circle")
      arr.width=arr.length
    rad <- 0.1                # radius of circle
    mid <- c(-rad,0)
    x<- seq(0,2*pi,length.out=15*npoint)
    rr <- 6.22*cbind(mid[1]+rad*sin(x),mid[2]+rad*cos(x))
  }

  if(arr.adj == 0.5)
    rr[,1] <- rr[,1]-min(rr[,1])/2
  if(arr.adj == 0)
    rr[,1] <- rr[,1]-min(rr[,1])

  user <- par("usr")
  pcm  <- par("pin")*2.54

  sy<- (user[4]-user[3])/pcm[2]
  sx<- (user[2]-user[1])/pcm[1]
  nr <- max(length(x0),length(y0),length(angle),
            length(arr.length),length(arr.width),
            length(lcol),length(lty),length(arr.col))
  if (nr>1) {
    x0         <- rep(x0        ,length.out=nr)
    y0         <- rep(y0        ,length.out=nr)
    angle      <- rep(angle     ,length.out=nr)
    arr.length <- rep(arr.length,length.out=nr)
    arr.width  <- rep(arr.width,length.out=nr)
    lcol       <- rep(lcol      ,length.out=nr)
    lty        <- rep(lty       ,length.out=nr)
    arr.col    <- rep(arr.col   ,length.out=nr)
  }
  RR<-rr
  for (i in 1:nr) {
  ## rotation around midpoint
    dx <- rr[,1]*arr.length [i]
    dy <- rr[,2]*arr.width  [i]

    angpi <- angle[i] / 180 *pi
    cosa  <-cos(angpi)
    sina  <-sin(angpi)

    RR[,1]<-  cosa*dx-sina*dy
    RR[,2]<-  sina*dx+cosa*dy

## rescaling and transposing
    RR[,1]<- x0[i] +RR[,1]*sx
    RR[,2]<- y0[i] +RR[,2]*sy

## drawing...
    polygon(RR,col=arr.col[i],border=lcol[i],lty=lty[i],lwd=arr.lwd)
  }
}
