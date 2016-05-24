`fun.data.fit.qs` <-
function (data, rs.leap = 3, fmkl.leap = 3, rs.init = c(-1.5, 
    1.5), fmkl.init = c(-0.25, 1.5), FUN = "runif.sobol",trial.n=100,len=1000,type=7,no=10000) 
{
    RPRS <- fun.RPRS.qs(data = data, rs.init = rs.init, leap = rs.leap, 
        FUN = FUN,trial.n=trial.n,len=len,type=type,no=no)
    RMFMKL <- fun.RMFMKL.qs(data = data, fmkl.init = fmkl.init, 
        leap = fmkl.leap, FUN = FUN,trial.n=trial.n,len=len,type=type,no=no)
    STAR <- starship(data)$lambda
    cbind(RPRS, RMFMKL)
}

