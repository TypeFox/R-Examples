`fun.RMFMKL.mm` <-
function (data, fmkl.init = c(-0.25, 1.5), leap = 3, FUN = "runif.sobol",no=10000) 
{
    RMFMKL <- fun.fit.gl.v4(a=fmkl.init[1], b=fmkl.init[2], data=data, fun=fun.auto.mm.fmkl,no=no,
leap = leap, FUN = FUN)$unique.optim.result
    RMFMKL <- fun.fit.gl.v4a(RMFMKL[1], RMFMKL[2], RMFMKL[3], 
        RMFMKL[4], data, "fmkl")
    return(RMFMKL)
}

