`fun.data.fit.lm` <-
function (data, rs.leap = 3, fmkl.leap = 3, rs.init = c(-1.5, 
    1.5), fmkl.init = c(-0.25, 1.5), FUN = "runif.sobol",no=10000) 
{
    RPRS <- fun.RPRS.lm(data = data, rs.init = rs.init, leap = rs.leap, 
        FUN = FUN,no=no)
    RMFMKL <- fun.RMFMKL.lm(data = data, fmkl.init = fmkl.init, 
        leap = fmkl.leap, FUN = FUN,no=no)
    cbind(RPRS, RMFMKL)
}

