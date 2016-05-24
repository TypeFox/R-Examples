`fun.RPRS.qs` <-
function (data, rs.init = c(-1.5, 1.5), leap = 3, FUN = "runif.sobol",trial.n=100,len=1000,type=7,no=10000) 
{
    RPRS <- fun.fit.gl.v6(a=rs.init[1], b=rs.init[2], data=data, fun=fun.auto.perc.rs, no=no,
        leap = leap, FUN = FUN,trial.n=trial.n,len=len,type=type)$unique.optim.result
    RPRS <- fun.fit.gl.v6a(RPRS[1], RPRS[2], RPRS[3], RPRS[4], 
        data, "rs",len=len,type=type)
    return(RPRS)
}

