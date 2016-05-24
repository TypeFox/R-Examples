DensCircRad<-function (x, t, bw) {
    exptx<-exp(cos(outer(t,x,"-"))*bw)
    y<-rowMeans(exptx)/(2*pi*besselI(bw,0))
    return(y)
}