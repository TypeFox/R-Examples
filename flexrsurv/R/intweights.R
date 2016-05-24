intweights_CAV_SIM<-function(nl){
  # wehights for cavalieri simpson rule (nl even)
  # weights are (1 4 2 4 2 4 ... 4 2 4 2 4 1 )
  w<-rep(2L/3L, nl+1)
  w[2*(1:(nl/2))]<-4L/3L
  w[1]<-w[nl+1]<-1L/3L
  w
}
  

intweights_SIM_3_8<-function(nl){
  # wehights for 3/8 simpson rule (nl = 3 * i)
  # weights are (1 3 3 2 3 3 2 3 3  ... 3 3 2 3 3 2 3 3 1 )
  w<-rep(9L/8L, nl+1)
  w[1+3*(1:(nl/3-1))]<-3L/4L
  w[1]<-w[nl+1]<-3L/8L
  w 
}


intweights_BOOLE<-function(nl){
  # wehights for BOOLE rule (nl = 4 * i)
  # weights are (7 32 12 32 14 32 12 32 14 ... 14 32 12 32 14 32 12 32 7 )
  w<-rep(64L/45L, nl+1)
  w[4*1:(nl/4)-1]<-8L/15L
  w[4*1:(nl/4)+1]<-28L/45L
  w[1]<-w[nl+1]<-14L/45L
  w
}
