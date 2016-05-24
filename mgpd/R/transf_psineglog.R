transf_psineglog <-
function 
(mar1 = c(0, 1, 0.1), mar2 = c(0, 1, 0.1), dep = 2, asy = 0, p = 3, compare=2,...)
{
    A1    = expression(1 - (x^(-alpha) + (1 - x)^(-alpha))^(-1/alpha))
    fi1   = expression(c * t^a * (1 - t)^a + t)
    
    d1A1  = D(A1,"x")
    d2A1  = D(d1A1,"x")
    A     = function(x, alpha) eval({x<-x; alpha<-alpha; A1})
    d1A   = function(x, alpha) eval({x<-x; alpha<-alpha; d1A1})
    d2A   = function(x, alpha) eval({x<-x; alpha<-alpha; d2A1})
    d1fi1 = D(fi1,"t")
    d2fi1 = D(d1fi1,"t")
    fi    = function(t, c, a) eval({t<-t; c<-c; a<-a; fi1})
    d1fi  = function(t, c, a) eval({t<-t; c<-c; a<-a; d1fi1})
    d2fi  = function(t, c, a) eval({t<-t; c<-c; a<-a; d2fi1})

    Afi   = function(t, alpha, c, a) A(fi(t, c, a), alpha)
    d1Afi = function(t, alpha, c, a) d1A(fi(t, c, a), alpha) * d1fi(t, c, a)
    d2Afi = function(t, alpha, c, a) d2A(fi(t, c, a), alpha) * (d1fi(t, c, a))^2 + d1A(fi(t, c, a), alpha) * d2fi(t,c, a)
    mu    = function(x, y, alpha, c, a) (1/x + 1/y) * Afi(x/(x + y), alpha, c, a)

    param = as.numeric(c(mar1, mar2, dep, asy, p))
    mux   = param[1]; muy   = param[4]
    sigx  = param[2]; sigy  = param[5]
    gamx  = param[3]; gamy  = param[6]
    alpha = param[7]
    asy   = param[8]; p     = param[9]

    hxy   = NULL
    error = FALSE
    xx    = seq(0, 1, 0.01)
    
    #if(min(d1Afi(xx,alpha,asy,p),na.rm=TRUE)<(-1) | max(d1Afi(xx,alpha,asy,p),na.rm=TRUE)>(+1)) error=TRUE
    if(min(d1Afi(xx,alpha,asy,p),na.rm=TRUE)<(-1) | max(d1Afi(xx,alpha,asy,p),na.rm=TRUE)>(+1)) error=TRUE
    if(min(d2Afi(xx,alpha,asy,p),na.rm=TRUE)<0 ) error=TRUE
    if(sigx<0 | sigy<0 | alpha>5 | alpha < 1.1 ) error=TRUE
    #if (min(d2Afi(xx, alpha, asy, p), na.rm = TRUE) < 0 | sigx < 0 | sigy < 0 | alpha < 1.0001 ) error = TRUE
    
    par(mfrow=c(1,2))
    plot(xx,fi(xx,asy,p)-xx,t="l",main=expression(c * x^a * (1 - x)^a),xlab="x",ylab="")
    abline(h=0,lty=3)
    plot(xx,d2Afi(xx,alpha,asy,p),t="l",main="Spectral density",xlab="x",ylab="",ylim=c(0,4))
    lines(xx,d2A(xx,compare),lty=3)
    spdens=cbind(xx,fi(xx,asy,p)-xx,d2Afi(xx,alpha,asy,p),d2A(xx,compare))
    spdens
}
