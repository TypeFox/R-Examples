`fun.RPRS.hs` <-
function (data, default = "Y", rs.init = c(-1.5, 1.5), no.c.rs = 50, 
    leap = 3, FUN = "runif.sobol",no=10000) 
{
    if (default == "Y") {
        no.c.rs <- fun.nclass.e(data)
    }
    RPRS <- fun.fit.gl.v2a(a=rs.init[1], b=rs.init[2], data=data, fun=fun.auto.perc.rs, no=no,
        nclass = no.c.rs, leap = leap, FUN = FUN)$unique.optim.result
    RPRS <- fun.fit.gl.v2b(RPRS[1], RPRS[2], RPRS[3], RPRS[4], 
        data, fun.auto.perc.rs, nclass = no.c.rs)
    return(RPRS)
}

