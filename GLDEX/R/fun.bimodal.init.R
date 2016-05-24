`fun.bimodal.init` <-
function (data1, data2, rs.leap1, fmkl.leap1, rs.init1, fmkl.init1, 
    rs.leap2, fmkl.leap2, rs.init2, fmkl.init2, fun1 = "runif.sobol", 
    fun2 = "runif.sobol",no=10000) 
{
    prop <- length(data1)/(length(data1) + length(data2))
    first.fit <- fun.data.fit.ml(data1, rs.leap = rs.leap1, fmkl.leap = fmkl.leap1, 
        rs.init = rs.init1, fmkl.init = fmkl.init1, FUN = fun1,no=no)
    second.fit <- fun.data.fit.ml(data2, rs.leap = rs.leap2, 
        fmkl.leap = fmkl.leap2, rs.init = rs.init2, fmkl.init = fmkl.init2, 
        FUN = fun2,no=no)
    return(list("prop"=prop, "first.fit"=first.fit, "second.fit"=second.fit))
}

