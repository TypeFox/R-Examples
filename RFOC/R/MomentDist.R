MomentDist<-function(E1, E2)
  {
    ###   given two right-hand Moment tensors,
    ###   get distance between them
    ###  according to Tape and Tape, 2012
BU = t(E1$vectors) %*% E2$vectors

ww = 0.5*sqrt(1+BU[1,1]+BU[2,2]+BU[3,3])
qw = c(ww, (BU[3,2] - BU[2,3])/(4*ww),(BU[1,3] - BU[3,1])/(4*ww),(BU[2,1] - BU[1,2])/(4*ww))

zw = qw[which.max(abs(qw))]

xi0 = 2*acos(abs(zw) )*180/pi
 
return(xi0)
###  angle in degrees
}



