plotor <-
function(x, y, ...){
    x<-as.numeric(x)
    y<-as.numeric(y)
    m<-length(x)
    n<-length(y)
    z<-sort(c(x,y))
    tl<-seq(max(min(x),min(y))*1.1, min(max(x),max(y))*.9, len=min(m,n))
    Fx<-ecdf(x)
    Gy<-ecdf(y)
    xt<-"t"
    yt<-"Emipirical Odds Ratio"
    eor.n<-((1-Fx(tl))*Gy(tl))
    eor.d<-((1-Gy(tl))*Fx(tl))
    eor<-((1-Fx(tl))*Gy(tl))/((1-Gy(tl))*Fx(tl))
    plot(tl[eor.n>0 & eor.d>0],eor[eor.n>0 & eor.d>0], xlab=xt, ylab=yt, ...)
}

