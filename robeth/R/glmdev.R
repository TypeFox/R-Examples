glmdev <-
function(y, ni, ci, wa, vtheta, offset=0, icase = .dFvGet()$ics)
{
        if(missing(y)) messagena("y")
        n <- length(y)
        if(missing(ni)) ni <- integer(n) 
        if(missing(ci)) messagena("ci")
        if(missing(wa)) messagena("wa") 
        if(missing(vtheta)) messagena("vtheta") 
        if (length(offset)==1) offset <- rep(0,n)
        dev <- single(1)
        thetas <- single(n)
        li  <- single(n)      
        sc  <- single(n)  
#       sink("GLMini.tmp")
        f.res <- .Fortran("glmdev",
                y = to.single(y),
                ni = to.integer(ni),
                ci = to.single(ci),
                wa = to.single(wa),
                vtheta = to.single(vtheta),
                oi = to.single(offset),
                n = to.integer(n),
                icase = to.integer(icase),
                dev = to.single(dev),
                thetas = to.single(thetas),
                li = to.single(li),
                sc = to.single(sc))
#       sink() 
        sc <- f.res$sc 
#       list(dev = f.res$dev, thetas = f.res$thetas, cis = f.res$cis,
#               was = f.res$was, li = f.res$li, sc = sc)
        list(dev = f.res$dev, thetas = f.res$thetas, li = f.res$li, sc = sc)
}
