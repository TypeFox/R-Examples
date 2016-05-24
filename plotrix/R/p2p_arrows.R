p2p_arrows<-function(x1,y1,x2,y2,space=0.05,col=par("fg"),...) {
 xspace<-(x2-x1)*space
 yspace<-(y2-y1)*space
 arrows(x1+xspace,y1+yspace,x2-xspace,y2-yspace,...)
}
