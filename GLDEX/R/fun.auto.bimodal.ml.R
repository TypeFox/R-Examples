"fun.auto.bimodal.ml"<-
function (data, per.of.mix = 0.01, clustering.m = clara, init1.sel = "rprs", 
    init2.sel = "rprs", init1=c(-1.5, 1.5), init2=c(-1.5, 1.5), leap1=3, leap2=3, fun1 = "runif.sobol", 
    fun2 = "runif.sobol", no = 10000,max.it=5000,optim.further="Y") 
{
    data.mod <- fun.class.regime.bi(data, per.of.mix, clustering.m)
    data1 <- data.mod$data.a
    data2 <- data.mod$data.b
    prop <- length(data1)/(length(data1) + length(data2))
    if (init1.sel == "rprs") {
        selc1 <- "rs"
        param1 <- "rs"
        first.fit <- fun.RPRS.ml(data = data1, rs.init = init1, 
            leap = leap1, FUN = fun1, no = no)
    }
    if (init1.sel == "rmfmkl") {
        selc1 <- "fmkl"
        param1 <- "fmkl"
        first.fit <- fun.RMFMKL.ml(data = data1, fmkl.init = init1, 
            leap = leap1, FUN = fun1, no = no)
    }
    if (init1.sel == "star") {
        selc1 <- "fmkl"
        param1 <- "fmkl"
        first.fit <- starship(data1)$lambda
    }
    if (init2.sel == "rprs") {
        selc2 <- "rs"
        param2 <- "rs"
        second.fit <- fun.RPRS.ml(data = data2, rs.init = init2, 
            leap = leap2, FUN = fun2, no = no)
    }
    if (init2.sel == "rmfmkl") {
        selc2 <- "fmkl"
        param2 <- "fmkl"
        second.fit <- fun.RMFMKL.ml(data = data2, fmkl.init = init2, 
            leap = leap2, FUN = fun2, no = no)
    }
    if (init2.sel == "star") {
        selc2 <- "fmkl"
        param2 <- "fmkl"
        second.fit <- starship(data2)$lambda
    }
    result <- optim(c(first.fit, second.fit, prop), optim.fun5, 
        data = data, param1 = param1, param2 = param2, control = list(maxit = max.it))
    result <- optim(c(result$par), optim.fun5, data = data, param1 = param1, 
        param2 = param2, control = list(maxit = max.it))

if(optim.further=="Y"){

result<-optim(c(result$par), optim.fun.bi.final, data = data, param1 = param1, 
param2 = param2, control = list(maxit = max.it))

}

    return(result)
}
