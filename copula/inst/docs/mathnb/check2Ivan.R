library(copula, lib.loc="../../../../copula.Rcheck")
source("../../../R/claytonExpr.R")
source("../../../R/gumbelExpr.R")
source("../../../R/frankExpr.R")
source("../../../R/plackettExpr.R")

source("../../../R/E.R")

source("../../../R/derCdfPdf.R")


set.seed(1234)

#### clayton
cop <- claytonCopula(2, dim=5)
u <- rCopula(10, cop)

dCdtheta(cop, u)
dCdu(cop, u)
derPdfWrtParams(cop, u)
derPdfWrtArgs(cop, u)

#### gumbel
cop <- gumbelCopula(2, dim=5)
u <- rCopula(10, cop)

dCdtheta(cop, u)
dCdu(cop, u)
derPdfWrtParams(cop, u)
derPdfWrtArgs(cop, u)


#### frank
cop <- frankCopula(2, dim=5)
u <- rCopula(10, cop)

dCdtheta(cop, u)
dCdu(cop, u)
derPdfWrtParams(cop, u)
derPdfWrtArgs(cop, u)

#### plackett
cop <- plackettCopula(2) # dim = 2 only
u <- rCopula(10, cop)

dCdtheta(cop, u)
dCdu(cop, u)
derPdfWrtParams(cop, u)
derPdfWrtArgs(cop, u)
