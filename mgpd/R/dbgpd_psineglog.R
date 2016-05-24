dbgpd_psineglog <-
function 
(x, y,
mar1 = c(0, 1, 0.1), mar2 = c(0, 1, 0.1),
dep = 2, asy = 0, p = 3, ...)
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
    
    d2Axx                   = d2Afi(xx,alpha,asy,p)
    d2Axx[d2Axx==-Inf]      = NA
    if(min(d2Axx,na.rm=TRUE)<0 ) error=TRUE
    if(sigx<0 | sigy<0 | alpha<0 | alpha > 6 ) error=TRUE
    
    if (!error)
    {
    tx    = (1 + gamx * (x - mux)/sigx)^(1/gamx)
    ty    = (1 + gamy * (y - muy)/sigy)^(1/gamy)
    tx0   = (1 + gamx * (-mux)/sigx)^(1/gamx)
    ty0   = (1 + gamy * (-muy)/sigy)^(1/gamy)
    dtx   = (1/sigx) * pmax((1 + gamx * (x - mux)/sigx), 0)^(1/gamx -1)
    dty   = (1/sigy) * pmax((1 + gamy * (y - muy)/sigy), 0)^(1/gamy -1)
    c0    = -mu(tx0, ty0, alpha, asy, p)
    mu1   = tx/(tx + ty)
    dxmu1 = ty/(tx + ty)^2
    dymu1 = (-tx)/(tx + ty)^2
    dxdymu1 = (tx - ty)/(tx + ty)^3
    dxdymu = (-1) * d1Afi(mu1, alpha, asy, p) * (dymu1/tx^2 + dxmu1/ty^2) + (1/tx + 1/ty) * (d2Afi(mu1, alpha,asy, p) * dxmu1 * dymu1 + d1Afi(mu1, alpha, asy,p) * dxdymu1)
    hxy   = 1/c0 * dxdymu * dtx * dty
    hxy   = as.numeric(hxy * (1 - ((x < 0) * (y < 0))))
    hxy
    }else stop("invalid parameter(s)")
    hxy
}
