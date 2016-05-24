`fun.data.fit.ml` <-
function (data, rs.leap = 3, fmkl.leap = 3, rs.init = c(-1.5, 
    1.5), fmkl.init = c(-0.25, 1.5), FUN = "runif.sobol",no=10000) 
{
    RPRS <- fun.RPRS.ml(data = data, rs.init = rs.init, leap = rs.leap, 
        FUN = FUN,no=no)
    RMFMKL <- fun.RMFMKL.ml(data = data, fmkl.init = fmkl.init, 
        leap = fmkl.leap, FUN = FUN,no=no)
    STAR <- starship(data)$lambda
    cbind(RPRS, RMFMKL, STAR)
}

