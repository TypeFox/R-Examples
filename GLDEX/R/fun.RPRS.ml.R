`fun.RPRS.ml` <-
function (data, rs.init = c(-1.5, 1.5), leap = 3, FUN = "runif.sobol",no=10000) 
{
    RPRS <- fun.fit.gl.v3(a=rs.init[1], b=rs.init[2], data=data, fun=fun.auto.perc.rs, no=no,
        leap = leap, FUN = FUN)$unique.optim.result
    RPRS <- fun.fit.gl.v3a(RPRS[1], RPRS[2], RPRS[3], RPRS[4], 
        data, "rs")
    return(RPRS)
}

