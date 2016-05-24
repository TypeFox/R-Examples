`fun.RMFMKL.qs` <-
function (data, fmkl.init = c(-0.25, 1.5), leap = 3, FUN = "runif.sobol",trial.n=100,len=1000,type=7,no=10000) 
{
    RMFMKL <- fun.fit.gl.v6(a=fmkl.init[1], b=fmkl.init[2], data=data, 
        fun=fun.auto.mm.fmkl, no=no, leap = leap, FUN = FUN,trial.n=trial.n,len=len,type=type)$unique.optim.result
    RMFMKL <- fun.fit.gl.v6a(RMFMKL[1], RMFMKL[2], RMFMKL[3], 
        RMFMKL[4], data, "fmkl",len=len,type=type)
    return(RMFMKL)
}

