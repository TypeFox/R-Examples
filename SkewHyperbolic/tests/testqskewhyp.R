require(SkewHyperbolic)
## source("../R/skewhypCalcRange.R")
## source("../R/dskewhyp.R")
## source("../data/skewhypParam.R")
## library(RUnit)

options(digits=20)
param <- c(0,1,0,10)
q <- c(-Inf,-1,0,1,Inf)
pskewhyp(q, param = param)
pskewhyp(q, param = param, lower.tail = FALSE)
pskewhyp(q, param = param, valueOnly = FALSE)
pskewhyp(q, param = param, valueOnly = FALSE, intTol = 10^(-12))

x <- rskewhyp(1, param = param)
x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-10)),
            param = param) - x
qskewhyp(pskewhyp(x, param = param),
            param = param, uniTol = 10^(-10)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-8)),
            param = param, uniTol = 10^(-8)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-10)),
            param = param, uniTol = 10^(-10)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-12)),
            param = param, uniTol = 10^(-12)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-10)),
         param = param, method = "integrate") - x
qskewhyp(pskewhyp(x, param = param),
            param = param, uniTol = 10^(-10)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-8)),
         param = param, uniTol = 10^(-8), method = "integrate") - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-10)),
         param = param, uniTol = 10^(-10), method = "integrate") - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-12)),
         param = param, uniTol = 10^(-12), method = "integrate") - x
qskewhyp(pskewhyp(10, param = param, intTol = 10^(-12)),
         param = param, uniTol = 10^(-12), method = "integrate") - 10
qskewhyp(pskewhyp(-10, param = param, intTol = 10^(-12)),
         param = param, uniTol = 10^(-12), method = "integrate") + 10


param <- c(0,1,10,20)
q <- c(-Inf,-1,0,1,Inf)
pskewhyp(q, param = param)
pskewhyp(q, param = param, lower.tail = FALSE)
pskewhyp(q, param = param, valueOnly = FALSE)
pskewhyp(q, param = param, valueOnly = FALSE, intTol = 10^(-12))

x <- rskewhyp(1, param = param)
x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-10)),
            param = param) - x
qskewhyp(pskewhyp(x, param = param),
            param = param, uniTol = 10^(-10)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-8)),
            param = param, uniTol = 10^(-8)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-10)),
            param = param, uniTol = 10^(-10)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-12)),
            param = param, uniTol = 10^(-12)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-10)),
         param = param, method = "integrate") - x
qskewhyp(pskewhyp(x, param = param),
            param = param, uniTol = 10^(-10)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-8)),
         param = param, uniTol = 10^(-8), method = "integrate") - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-10)),
         param = param, uniTol = 10^(-10), method = "integrate") - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-12)),
         param = param, uniTol = 10^(-12), method = "integrate") - x
qskewhyp(pskewhyp(10, param = param, intTol = 10^(-12)),
         param = param, uniTol = 10^(-12), method = "integrate") - 10
qskewhyp(pskewhyp(-10, param = param, intTol = 10^(-12)),
         param = param, uniTol = 10^(-12), method = "integrate") + 10

param <- c(0,1,1,1)
q <- c(-Inf,-1,0,1,Inf)
pskewhyp(q, param = param)
pskewhyp(q, param = param, lower.tail = FALSE)
pskewhyp(q, param = param, valueOnly = FALSE)
pskewhyp(q, param = param, valueOnly = FALSE, intTol = 10^(-12))


param <- c(0,1,-10,5)
q <- c(-Inf,-1,0,1,Inf)
pskewhyp(q, param = param)
pskewhyp(q, param = param, lower.tail = FALSE)
pskewhyp(q, param = param, valueOnly = FALSE)
pskewhyp(q, param = param, valueOnly = FALSE, intTol = 10^(-12))



param <- c(0,1,5,5)
q <- c(-Inf,-1,0,1,Inf)
pskewhyp(q, param = param)
pskewhyp(q, param = param, lower.tail = FALSE)
pskewhyp(q, param = param, valueOnly = FALSE)
pskewhyp(q, param = param, valueOnly = FALSE, intTol = 10^(-12))

x <- rskewhyp(1, param = param)
x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-10)),
            param = param) - x
qskewhyp(pskewhyp(x, param = param),
            param = param, uniTol = 10^(-10)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-8)),
            param = param, uniTol = 10^(-8)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-10)),
            param = param, uniTol = 10^(-10)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-12)),
            param = param, uniTol = 10^(-12)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-10)),
         param = param, method = "integrate") - x
qskewhyp(pskewhyp(x, param = param),
            param = param, uniTol = 10^(-10)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-8)),
         param = param, uniTol = 10^(-8), method = "integrate") - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-10)),
         param = param, uniTol = 10^(-10), method = "integrate") - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-12)),
         param = param, uniTol = 10^(-12), method = "integrate") - x
qskewhyp(pskewhyp(10, param = param, intTol = 10^(-12)),
         param = param, uniTol = 10^(-12), method = "integrate") - 10
qskewhyp(pskewhyp(-10, param = param, intTol = 10^(-12)),
         param = param, uniTol = 10^(-12), method = "integrate") + 10

param <- c(1,2,20,10)
q <- c(-Inf,-1,0,1,Inf)
pskewhyp(q, param = param)
pskewhyp(q, param = param, lower.tail = FALSE)
pskewhyp(q, param = param, valueOnly = FALSE)
pskewhyp(q, param = param, valueOnly = FALSE, intTol = 10^(-12))

x <- rskewhyp(1, param = param)
x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-10)),
            param = param) - x
qskewhyp(pskewhyp(x, param = param),
            param = param, uniTol = 10^(-10)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-8)),
            param = param, uniTol = 10^(-8)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-10)),
            param = param, uniTol = 10^(-10)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-12)),
            param = param, uniTol = 10^(-12)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-10)),
         param = param, method = "integrate") - x
qskewhyp(pskewhyp(x, param = param),
            param = param, uniTol = 10^(-10)) - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-8)),
         param = param, uniTol = 10^(-8), method = "integrate") - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-10)),
         param = param, uniTol = 10^(-10), method = "integrate") - x
qskewhyp(pskewhyp(x, param = param, intTol = 10^(-12)),
         param = param, uniTol = 10^(-12), method = "integrate") - x
qskewhyp(pskewhyp(10, param = param, intTol = 10^(-12)),
         param = param, uniTol = 10^(-12), method = "integrate") - 10
qskewhyp(pskewhyp(-10, param = param, intTol = 10^(-12)),
         param = param, uniTol = 10^(-12), method = "integrate") + 10
