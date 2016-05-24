`fun.data.fit.hs.nw` <-
function (data, rs.default = "Y", fmkl.default = "Y", rs.leap = 3, 
    fmkl.leap = 3, rs.init = c(-1.5, 1.5), fmkl.init = c(-0.25, 
        1.5), no.c.rs = 50, no.c.fmkl = 50, FUN = "runif.sobol",no=10000) 
{
    RPRS.nw <- fun.RPRS.hs.nw(data = data, default = rs.default, 
        rs.init = rs.init, no.c.rs = no.c.rs, leap = rs.leap, 
        FUN = FUN,no=no)
    RMFMKL.nw <- fun.RMFMKL.hs.nw(data = data, default = fmkl.default, 
        fmkl.init = fmkl.init, no.c.fmkl = no.c.fmkl, leap = fmkl.leap, 
        FUN = FUN,no=no)
    cbind(RPRS.nw, RMFMKL.nw)
}

