disturber <- function(tt,y,s,ps, s.outlier.fraction=0, interval) {
    if(interval & (max(tt)-3*ps)<0) stop("Light curve is too short to allow interval=TRUE")
    n.points<- length(tt)
    if(s.outlier.fraction>0) {
        hnr <- sample((1:n.points), round(n.points*s.outlier.fraction), replace=FALSE)
        s[hnr]<- 0.5*min(s)
    }
    if(interval){
        peak.start <- sample(tt[which(tt< (max(tt)-3*ps))],1)
        peak.hnr <- which(tt>=peak.start & tt<= peak.start+3*ps)
        y[peak.hnr] <- 6*quantile(y,0.9)* dnorm(tt[peak.hnr], mean=peak.start+1.5*ps, sd=ps)/dnorm(0, sd=ps)
    }
    return(list(y=y, s=s))
}
