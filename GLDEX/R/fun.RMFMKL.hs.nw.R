`fun.RMFMKL.hs.nw` <-
function (data, default = "Y", fmkl.init = c(-0.25, 1.5), no.c.fmkl = 50, 
    leap = 3, FUN = "runif.sobol",no=10000) 
{
    if (default == "Y") {
        no.c.fmkl <- fun.nclass.e(data)
    }
    RMFMKL <- fun.fit.gl.v2a.nw(a=fmkl.init[1], b=fmkl.init[2], data=data, 
        fun=fun.auto.mm.fmkl, no=no,nclass = no.c.fmkl, leap = leap, FUN = FUN)$unique.optim.result
    RMFMKL <- fun.fit.gl.v2b.nw(RMFMKL[1], RMFMKL[2], RMFMKL[3], 
        RMFMKL[4], data, fun.auto.mm.fmkl, nclass = no.c.fmkl)
    return(RMFMKL)
}

