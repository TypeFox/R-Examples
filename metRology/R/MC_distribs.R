
#Distribution functions
.urnorm<-function(B, x, u) {
        rnorm(B, x, u)

}

.urunif<-function(B, x, u) {
        runif(B, x-u*sqrt(3), x+u*sqrt(3))
}

#Triangular
dtri<-function(x, min=-sqrt(6), max=sqrt(6), mode=(min+max)/2, log=FALSE) {
        lo <- pmin(min,max)
        hi <- pmax(min,max)
        if(any((mode < lo) | (mode > hi) )) warning("Mode outside interval", call.=TRUE)
        
        d<-ifelse(x <= mode, 2*(x-min) / ((hi-lo)*(mode-lo)), 2*(hi-x) / ((hi-lo)*(hi-mode)))
        d[x < lo | x > hi] <- 0
        d[(mode < lo) | (mode > hi)] <- NA
        
        if(log) return(log(d)) else return(d)
}

ptri<-function(q, min=-sqrt(6), max=sqrt(6), mode=(min+max)/2, lower.tail = TRUE, log.p=FALSE) {
        lo <- pmin(min,max)
        hi <- pmax(min,max)
        if(any((mode < lo) | (mode > hi) )) warning("Mode outside interval", call.=TRUE)
        
        p<-ifelse(q <= mode, (q-lo)^2 / ((mode-lo)*(hi-lo)), 1-(hi-q)^2 / ((hi-mode)*(hi-lo)))
        p[q < lo] <- 0
        p[ q > hi] <- 1
        p[(mode < lo) | (mode > hi)] <- NA
        
        if(!lower.tail) p <- 1-p
        if(log.p) return(log(p)) else return(p)
}

qtri<-function(p, min=-sqrt(6), max=sqrt(6), mode=(min+max)/2, lower.tail = TRUE, log.p=FALSE) {
        if(log.p) p <- exp(p)
        if(!lower.tail) p <- 1-p
        lo <- pmin(min,max)
        hi <- pmax(min,max)
        if(any((mode < lo) | (mode > hi) )) warning("Mode outside interval", call.=TRUE)
        if(any((p < 0) | (p > 1) )) warning("p must be in the interval (0,1)", call.=TRUE)
        
        q<-ifelse(p <= (mode-lo)/(hi-lo), lo + sqrt(p*(hi-lo)*(mode-lo)), hi-sqrt((1-p)*(hi-lo)*(hi-mode)))
        q[p < 0] <- NA
        q[ p > 1] <- NA
        q[(mode < lo) | (mode > hi)] <- NA

        return(q)
        
}

rtri<-function(n, min=-sqrt(6), max=sqrt(6), mode=(min+max)/2) {
        p<-runif(n,min=0,max=1)
        return(qtri(p, min,max,mode))
}


#t
dt.scaled<-function(x, df, mean=0, sd=1, ncp, log = FALSE) {
        stats::dt((x-mean)/sd, df, ncp=ncp, log=FALSE)/sd
}
pt.scaled<-function(q, df, mean=0, sd=1, ncp, lower.tail = TRUE, log.p = FALSE) {
        stats::pt((q-mean)/sd, df, ncp=ncp, log.p=FALSE)
}
qt.scaled<-function(p, df, mean=0, sd=1, ncp, lower.tail = TRUE, log.p = FALSE) {
        mean+sd*stats::qt(p, df, ncp=ncp, log.p=FALSE)
}
rt.scaled<-function(n, df, mean=0, sd=1, ncp) {
        mean+sd*stats::rt(n, df, ncp=ncp)
}

