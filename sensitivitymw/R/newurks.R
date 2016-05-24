newurks <-
function (smat, m1 = 2, m2 = 2, m = 2) 
    {
        rg <- apply(smat, 1, max) - apply(smat, 1, min)
        rk <- rank(rg)
        urk <- multrnks(rk, m1 = m1, m2 = m2, m = m)
        for (j in 1:(dim(smat)[2])) {
            smat[, j] <- smat[, j] * urk
        }
        smat
    }

