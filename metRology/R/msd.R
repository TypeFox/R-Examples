#Calculates median scaled difference for observations x given sd's s
#if s is a scalar function, it is applied to x to obtain an estimate of s

msd<-function(x, s=mad , ...) {
        ss <- if(is.function(s)) {
                rep(s(x, ...), length(x))
        } else {
                if(length(s) == 1) {
                        rep(s, length(x))
                } else {
                        s
                }
        }
        
        m<-abs(outer(x,x,"-")) / outer(ss,ss,FUN=function(a,b) sqrt(a^2+b^2))
        diag(m) <- NA   #removes deviation from self
        
        return(apply(m, 1, median, na.rm=TRUE))
        
}

.pmsd.xnorm<-function(q, x, n, sd=1, scale=FALSE) {
        #P for the median of abs(x-X) (optionally /sqrt(2))
        #based on Mood, Graybill and Boes (1974) pp252ff
        #q can be a vector, n can't
        #NOTE: n is the number of values, not the number of differences.
        
        if(scale) sd <- sd/sqrt(2)
        
        pxnorm<-function(q,x,sd=1) ifelse(q>0, pnorm(x+q, 0, sd)-pnorm(x-q, 0, sd), 0) 
        
        Fy<-rep(0, length(q))
        
        n.med<-floor((n+1)/2)   #exact for odd samples, low for even
                                #Note that for n values there are n-1 differences,
                                #so an even-n case is an odd-median case
        
        n<-2*n.med              #Corrects n to the next higher even n
        
        ph<-pxnorm(q,x,sd)
        
        for(j in n.med:(n-1)) Fy <- Fy + choose(n-1,j) * (ph^j) * (1-ph)^(n-j-1)
        
        return(Fy)
        
}

pmsd<-function(q, n, sd=1, scale=TRUE) {
        #P for the median of abs(x-X)/sqrt(2) 
        #based on Mood, Graybill and Boes (1974) pp252ff
        #q can be a vector, n can't
        #NOTE:  In this implementation, n is the number of values, so msd 
        #       is the median of n-1 differences
        
        L <- max(length(q), length(n), length(sd))
        q<-rep(q, length.out=L)
        n<-rep(n, length.out=L)
        sd<-rep(sd, length.out=L)
        
        if(scale) sd <- sd/sqrt(2)

        f<-function(x, q, n, sd=1, scale=FALSE) {
                L<-length(x)
                q<-rep(q, length.out=L)
                n<-rep(n, length.out=L)
                sd <- rep(sd, length.out=L)
                rv<-rep(0, length(x))
                for(i in 1:length(rv)) rv[i]<-.pmsd.xnorm(q[i], x[i],  n[i], sd[i], scale)
                return(rv*dnorm(x,0,sd))
        }
        
        p <- rep(0, length(q))
        
        for(i in 1:length(q)) p[i]<-integrate(f, lower=-Inf, upper=Inf, rel.tol = .Machine$double.eps^0.75, 
                q=q[i], n=n[i], sd=sd[i], scale=FALSE)$value
                             #Note odd tolerance; pmsd and qmsd are quite inaccurate 
                             #(and qmsd even unstable) with default integrate() tolerance
        return(p)
        
}

qmsd<-function(p, n, sd=1, scale=TRUE) {
        #returns quantile(s) for msd using numerical root-finder
        #on pmsd.
        #p or n may be vectors
        #scale allows suppression of the scale parameter
        
        L <- max(length(q), length(n), length(sd))
        p<-rep(p, length.out=L)
        n<-rep(n, length.out=L)
        sd<-rep(sd, length.out=L)

        q<-rep(-1, length(p))
        
        q[p==0] <- 0
        q[p==1] <- +Inf

        froot<-function(q, p, n, sd=1, scale) pmsd(q, n, sd, scale)-p

        q.upper <-ifelse(p> 0.9999, 10,4) #saves a step
        
        for(i in 1:length(p)) {
                
                if(q[i]<0) { #Note explicit values above
                             #Also note odd tolerance; qmsd is quite inaccurate with default tolerance
                        q[i] <- uniroot(froot, interval=c(0,q.upper[i]), tol=.Machine$double.eps^0.75,
                                p=p[i], n=n[i], 
                                sd=sd[i], scale=scale)$root
                }
        }

        
        return(q)
}
