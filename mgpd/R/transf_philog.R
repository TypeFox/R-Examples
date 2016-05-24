transf_philog <-
function (
mar1 = c(0, 1, 0.1), mar2 = c(0, 1, 0.1),
dep = 2, asy = 0, p = 2, compare=2, ...)
{

fipoints=function(psi1,psi2)
{
points(0,0,cex=0.7)
points(1,0,cex=0.7)
points(1/psi2,0,cex=0.7)
points(1/psi2/2,psi1,cex=0.7)
points(1/2+1/psi2/2,-psi1,cex=0.7)
}
    A1    = expression((x^alpha + (1 - x)^alpha)^(1/alpha))
    fi1   = expression(t+((64*a*(c^5-c^4-2*c^3-2*c^2+5*c-2)*c^5)/((c+1)^2*(c-1)^2
*(-1+2*c)*(1+2*c^2-3*c)))*t^6+(-32*a*(5*c^6-2*c^5-13*c^4-12*c^3+15*c^2+6*c-6)*
c^4/((c+1)^2*(c-1)^2*(-1+2*c)*(1+2*c^2-3*c)))*t^5+(32*c^3*a*(4*c^7+3*c^6-14*c^5-
15*c^4+23*c^2-8*c-2)/((c+1)^2*(c-1)^2*(-1+2*c)*(1+2*c^2-3*c)))*t^4+(-(32*(c^7+4*
c^6-5*c^5-10*c^4-5*c^3+12*c^2+2*c-4))*a*c^3/((1+2*c^2-3*c)*(c+1)^2*
(4*c-1+2*c^3-5*c^2)))*t^3+((32*c^3*a*(c^6-3*c^4-c^2+4*c-2))/((1+2*c^2-3*c)*
(c+1)^2*(4*c-1+2*c^3-5*c^2)))*t^2)
#psi1=a ilyen magassagban megy at = ez kene legyen az asy
#psi2=c metszespont reciproka

    d1A1  = D(A1,"x")
    d2A1  = D(d1A1,"x")
    A     = function(x, alpha) eval({x<-x; alpha<-alpha; A1})
    d1A   = function(x, alpha) eval({x<-x; alpha<-alpha; d1A1})
    d2A   = function(x, alpha) eval({x<-x; alpha<-alpha; d2A1})
    d1fi1 = D(fi1,"t")
    d2fi1 = D(d1fi1,"t")
    fi    = function(t, a, c) eval({t<-t; c<-c; a<-a; fi1})
    d1fi  = function(t, a, c) eval({t<-t; c<-c; a<-a; d1fi1})
    d2fi  = function(t, a, c) eval({t<-t; c<-c; a<-a; d2fi1})

    Afi   = function(t, alpha, a, c) A(fi(t, a, c), alpha)
    d1Afi = function(t, alpha, a, c) d1A(fi(t, a, c), alpha) * d1fi(t, a, c)
    d2Afi = function(t, alpha, a, c) d2A(fi(t, a, c), alpha) * (d1fi(t, a, c))^2 + d1A(fi(t, a, c), alpha) * d2fi(t,a, c)
    mu    = function(x, y, alpha, a, c) (1/x + 1/y) * Afi(x/(x + y), alpha, a, c)

    param = as.numeric(c(mar1, mar2, dep, asy, p))
    mux   = param[1]; muy   = param[4]
    sigx  = param[2]; sigy  = param[5]
    gamx  = param[3]; gamy  = param[6]
    alpha = param[7]
    asy   = param[8]; p     = param[9]
    
    hxy   = NULL
    error = FALSE
    xx    = seq(0, 1, 0.005)
    
    par(mfrow=c(1,2))
    fixx  = fi(xx,asy,p)-xx
    plot(xx,fixx,t="l",ylim=c(min(fixx),max(fixx)),main="",xlab="x",ylab="")
    fipoints(asy,p)
    d2Axx = d2Afi(xx,alpha, asy, p)
    d2Axx[d2Axx>100]=NA
    plot(xx,d2Axx,t="l",ylim=c(0,4),main="Spectral density",xlab="x",ylab="")
    spdens=cbind(xx,fixx,d2Axx)
    abline(h=0,lty=2)
    d2Axx = d2Afi(xx,compare, 0, p)
    d2Axx[d2Axx>100]=NA
    lines(xx,d2Axx,lty=3)
    spdens
}
