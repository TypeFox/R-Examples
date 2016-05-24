pkgname <- "ESG"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('ESG')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Asset_PriceDistribution")
### * Asset_PriceDistribution

flush(stderr()); flush(stdout())

### Name: Asset_PriceDistribution
### Title: Asset_PriceDistribution method
### Aliases: Asset_PriceDistribution
###   Asset_PriceDistribution,Scenarios-method

### ** Examples

objScenario  <- new("Scenarios")
# Basic scenario's parameters setting
objScenario  <- setParamsBaseScenarios(objScenario, horizon = 10, nScenarios = 1000)
# Risk factors parameters setting
objScenario  <- setRiskParamsScenariosrt(objScenario, vol = .1, k = 2)
objScenario  <- setRiskParamsScenariosS(objScenario, vol = .1, k = 2, volStock = .2, stock0 = 100, rho=.5)
objScenario  <- setRiskParamsScenariosliqSpr(objScenario, eta=.05, liquiditySpread0=.01)
objScenario  <- setRiskParamsScenariosdefSpr(objScenario, volDefault=.2, defaultSpread0=.01, alpha=.1, beta=1)
# Forward and ZC rates setting
data(ZC)
objScenario  <- setForwardRates(objScenario, ZC, horizon=10)
objScenario  <- setZCRates(objScenario, ZC, horizon=10)
# Projection
objScenario  <- customPathsGeneration(objScenario, type="shortRate")
objScenario  <- customPathsGeneration(objScenario, type="stock")
objScenario  <- customPathsGeneration(objScenario, type="defaultSpread")
objScenario  <- customPathsGeneration(objScenario, type="liquiditySpread")
Asset_PriceDistribution(objScenario,type='ConvBond',t=0,T=10,nCoupons=1,couponsRate=0.03)



cleanEx()
nameEx("MartingaleTest")
### * MartingaleTest

flush(stderr()); flush(stdout())

### Name: MartingaleTest
### Title: MartingaleTest method
### Aliases: MartingaleTest MartingaleTest,Scenarios-method

### ** Examples

objScenario  <- new("Scenarios")
# Basic scenario's parameters setting
objScenario  <- setParamsBaseScenarios(objScenario, horizon = 10, nScenarios = 1000)
# Risk factors parameters setting
objScenario  <- setRiskParamsScenariosrt(objScenario, vol = .1, k = 2)
objScenario  <- setRiskParamsScenariosS(objScenario, vol = .1, k = 2, volStock = .2, stock0 = 100, rho=.5)
objScenario  <- setRiskParamsScenariosliqSpr(objScenario, eta=.05, liquiditySpread0=.01)
objScenario  <- setRiskParamsScenariosdefSpr(objScenario, volDefault=.2, defaultSpread0=.01, alpha=.1, beta=1)
# Forward and ZC rates setting
data(ZC)
objScenario  <- setForwardRates(objScenario, ZC, horizon=10)
objScenario  <- setZCRates(objScenario, ZC, horizon=10)
# Projection
objScenario  <- customPathsGeneration(objScenario, type="shortRate")
objScenario  <- customPathsGeneration(objScenario, type="stock")
objScenario  <- customPathsGeneration(objScenario, type="defaultSpread")
objScenario  <- customPathsGeneration(objScenario, type="liquiditySpread")
MartingaleTest(objScenario)



cleanEx()
nameEx("ZC")
### * ZC

flush(stderr()); flush(stdout())

### Name: ZC
### Title: ZC data
### Aliases: ZC
### Keywords: datasets

### ** Examples

data(ZC)



cleanEx()
nameEx("customPathsGeneration")
### * customPathsGeneration

flush(stderr()); flush(stdout())

### Name: customPathsGeneration
### Title: customPathsGeneration method
### Aliases: customPathsGeneration customPathsGeneration,Scenarios-method

### ** Examples

objScenario  <- new("Scenarios")
# Basic scenario's parameters setting
objScenario  <- setParamsBaseScenarios(objScenario, horizon = 10, nScenarios = 1000)
# Risk factors parameters setting
objScenario  <- setRiskParamsScenariosrt(objScenario, vol = .1, k = 2)
objScenario  <- setRiskParamsScenariosS(objScenario, vol = .1, k = 2, volStock = .2, stock0 = 100, rho=.5)
objScenario  <- setRiskParamsScenariosliqSpr(objScenario, eta=.05, liquiditySpread0=.01)
objScenario  <- setRiskParamsScenariosdefSpr(objScenario, volDefault=.2, defaultSpread0=.01, alpha=.1, beta=1)
# Forward and ZC rates setting
data(ZC)
objScenario  <- setForwardRates(objScenario, ZC, horizon=10)
objScenario  <- setZCRates(objScenario, ZC, horizon=10)
# Projection
objScenario  <- customPathsGeneration(objScenario, type="shortRate")
objScenario  <- customPathsGeneration(objScenario, type="stock")
objScenario  <- customPathsGeneration(objScenario, type="defaultSpread")
objScenario  <- customPathsGeneration(objScenario, type="liquiditySpread")



cleanEx()
nameEx("getForwardRates")
### * getForwardRates

flush(stderr()); flush(stdout())

### Name: getForwardRates
### Title: getForwardRates method
### Aliases: getForwardRates getForwardRates,Scenarios-method

### ** Examples

scenarios1 <- new("Scenarios")
scenarios1 <- setParamsBaseScenarios(scenarios1, horizon=5, nScenarios=10)
scenarios1 <- setRiskParamsScenarios(scenarios1, vol=.1, k=2,volStock=.2, volRealEstate=.15,volDefault=.2, alpha=.1,beta=1, eta=.05,rho=.5, stock0=100,realEstate0=50, liquiditySpread0=.01, defaultSpread0=.01)
data(ZC)
scenarios1 <- setForwardRates(scenarios1, ZC, horizon=5)
getForwardRates(scenarios1)



cleanEx()
nameEx("getParamsBaseScenarios")
### * getParamsBaseScenarios

flush(stderr()); flush(stdout())

### Name: getParamsBaseScenarios
### Title: getParamsBaseScenarios method
### Aliases: getParamsBaseScenarios getParamsBaseScenarios,Scenarios-method

### ** Examples

scenarios1 <- new("Scenarios")
scenarios1 <- setParamsBaseScenarios(scenarios1, horizon=5, nScenarios=10)
getParamsBaseScenarios(scenarios1)



cleanEx()
nameEx("getRiskParamsScenarios")
### * getRiskParamsScenarios

flush(stderr()); flush(stdout())

### Name: getRiskParamsScenarios
### Title: getRiskParamsScenarios method
### Aliases: getRiskParamsScenarios getRiskParamsScenarios,Scenarios-method

### ** Examples

scenarios1 <- new("Scenarios")
scenarios1 <- setParamsBaseScenarios(scenarios1, horizon=5, nScenarios=10)
scenarios1 <- setRiskParamsScenarios(scenarios1, vol=.1, k=2,volStock=.2, volRealEstate=.15,volDefault=.2, alpha=.1,beta=1, eta=.05,rho=.5, stock0=100,realEstate0=50, liquiditySpread0=.01, defaultSpread0=.01)
getRiskParamsScenarios(scenarios1)



cleanEx()
nameEx("getRiskParamsScenariosRE")
### * getRiskParamsScenariosRE

flush(stderr()); flush(stdout())

### Name: getRiskParamsScenariosRE
### Title: getRiskParamsScenariosRE method
### Aliases: getRiskParamsScenariosRE
###   getRiskParamsScenariosRE,Scenarios-method

### ** Examples

scenarios1 <- new("Scenarios")
scenarios1 <- setParamsBaseScenarios(scenarios1, horizon=5, nScenarios=10)
scenarios1 <- setRiskParamsScenarios(scenarios1, vol=.1, k=2,volStock=.2, volRealEstate=.15,volDefault=.2, alpha=.1,beta=1, eta=.05,rho=.5, stock0=100,realEstate0=50, liquiditySpread0=.01, defaultSpread0=.01)
getRiskParamsScenariosRE(scenarios1)



cleanEx()
nameEx("getRiskParamsScenariosS")
### * getRiskParamsScenariosS

flush(stderr()); flush(stdout())

### Name: getRiskParamsScenariosS
### Title: getRiskParamsScenariosS method
### Aliases: getRiskParamsScenariosS
###   getRiskParamsScenariosS,Scenarios-method

### ** Examples

scenarios1 <- new("Scenarios")
scenarios1 <- setParamsBaseScenarios(scenarios1, horizon=5, nScenarios=10)
scenarios1 <- setRiskParamsScenarios(scenarios1, vol=.1, k=2,volStock=.2, volRealEstate=.15,volDefault=.2, alpha=.1,beta=1, eta=.05,rho=.5, stock0=100,realEstate0=50, liquiditySpread0=.01, defaultSpread0=.01)
getRiskParamsScenariosS(scenarios1)



cleanEx()
nameEx("getRiskParamsScenariosdefSpr")
### * getRiskParamsScenariosdefSpr

flush(stderr()); flush(stdout())

### Name: getRiskParamsScenariosdefSpr
### Title: getRiskParamsScenariosdefSpr method
### Aliases: getRiskParamsScenariosdefSpr
###   getRiskParamsScenariosdefSpr,Scenarios-method

### ** Examples

scenarios1 <- new("Scenarios")
scenarios1 <- setParamsBaseScenarios(scenarios1, horizon=5, nScenarios=10)
scenarios1 <- setRiskParamsScenarios(scenarios1, vol=.1, k=2,volStock=.2, volRealEstate=.15,volDefault=.2, alpha=.1,beta=1, eta=.05,rho=.5, stock0=100,realEstate0=50, liquiditySpread0=.01, defaultSpread0=.01)
getRiskParamsScenariosdefSpr(scenarios1)



cleanEx()
nameEx("getRiskParamsScenariosliqSpr")
### * getRiskParamsScenariosliqSpr

flush(stderr()); flush(stdout())

### Name: getRiskParamsScenariosliqSpr
### Title: getRiskParamsScenariosliqSpr method
### Aliases: getRiskParamsScenariosliqSpr
###   getRiskParamsScenariosliqSpr,Scenarios-method

### ** Examples

scenarios1 <- new("Scenarios")
scenarios1 <- setParamsBaseScenarios(scenarios1, horizon=5, nScenarios=10)
scenarios1 <- setRiskParamsScenarios(scenarios1, vol=.1, k=2,volStock=.2, volRealEstate=.15,volDefault=.2, alpha=.1,beta=1, eta=.05,rho=.5, stock0=100,realEstate0=50, liquiditySpread0=.01, defaultSpread0=.01)
getRiskParamsScenariosliqSpr(scenarios1)



cleanEx()
nameEx("getRiskParamsScenariosrt")
### * getRiskParamsScenariosrt

flush(stderr()); flush(stdout())

### Name: getRiskParamsScenariosrt
### Title: getRiskParamsScenariosrt method
### Aliases: getRiskParamsScenariosrt
###   getRiskParamsScenariosrt,Scenarios-method

### ** Examples

scenarios1 <- new("Scenarios")
scenarios1 <- setParamsBaseScenarios(scenarios1, horizon=5, nScenarios=10)
scenarios1 <- setRiskParamsScenarios(scenarios1, vol=.1, k=2,volStock=.2, volRealEstate=.15,volDefault=.2, alpha=.1,beta=1, eta=.05,rho=.5, stock0=100,realEstate0=50, liquiditySpread0=.01, defaultSpread0=.01)
getRiskParamsScenariosrt(scenarios1)



cleanEx()
nameEx("getZCRates")
### * getZCRates

flush(stderr()); flush(stdout())

### Name: getZCRates
### Title: getZCRates method
### Aliases: getZCRates getZCRates,Scenarios-method

### ** Examples

scenarios1 <- new("Scenarios")
scenarios1 <- setParamsBaseScenarios(scenarios1, horizon=5, nScenarios=10)
scenarios1 <- setRiskParamsScenarios(scenarios1, vol=.1, k=2,volStock=.2, volRealEstate=.15,volDefault=.2, alpha=.1,beta=1, eta=.05,rho=.5, stock0=100,realEstate0=50, liquiditySpread0=.01, defaultSpread0=.01)
data(ZC)
scenarios1 <- setZCRates(scenarios1, ZC, horizon=5)
getZCRates(scenarios1)



cleanEx()
nameEx("rAllRisksFactors")
### * rAllRisksFactors

flush(stderr()); flush(stdout())

### Name: rAllRisksFactors
### Title: rAllRisksFactors
### Aliases: rAllRisksFactors

### ** Examples

data(ZC)
rAllRisksFactors(horizon=5, nScenarios=10, ZC, vol=.1, k=2, volStock=.2, stock0=100, rho=.5, volRealEstate=.15, realEstate0=50, eta=.05, liquiditySpread0=.01, defaultSpread0=.01, volDefault=.2, alpha=.1, beta=1)



cleanEx()
nameEx("rAssetDistribution")
### * rAssetDistribution

flush(stderr()); flush(stdout())

### Name: rAssetDistribution
### Title: rAssetDistribution
### Aliases: rAssetDistribution

### ** Examples

data(ZC)
rAssetDistribution(type="Zero-Coupon",t=2,T=3,vol=.1, k=2, ZC=ZC, nScenarios=100)
rAssetDistribution(type="Bond",t=3,T=35,nCoupons=20, couponsRate=0.3,vol=.1, k=2, ZC=ZC, nScenarios=10)
rAssetDistribution(type="CBond",t=5,T=35,nCoupons=5, couponsRate=0.3, omega=5,vol=.1, k=2, ZC=ZC, nScenarios=10,eta=.05, liquiditySpread0=.01, defaultSpread0=.01, volDefault=.2, alpha=.1, beta=1)
rAssetDistribution(type="EuroPut_Stock",5,25,Strike=98.5,vol=.1,k=2,ZC=ZC,volStock=.2, stock0=100, rho=.5,nScenarios=10)
rAssetDistribution(type="EuroCall_ZC",4,4.5,s=5, Strike=.985,vol=.1, k=2, ZC=ZC,nScenarios=10)
rAssetDistribution(type="EuroPut_ZC",4,4.5,s=5, Strike=.9385,vol=.1, k=2, ZC=ZC,nScenarios=10)



cleanEx()
nameEx("rDefaultSpread")
### * rDefaultSpread

flush(stderr()); flush(stdout())

### Name: rDefaultSpread
### Title: rDefaultSpread
### Aliases: rDefaultSpread

### ** Examples

rDefaultSpread(horizon=5, nScenarios=8, defaultSpread0=.01, volDefault=.2, alpha=.1, beta=1)



cleanEx()
nameEx("rLiquiditySpread")
### * rLiquiditySpread

flush(stderr()); flush(stdout())

### Name: rLiquiditySpread
### Title: rLiquiditySpread
### Aliases: rLiquiditySpread

### ** Examples

rLiquiditySpread(horizon=5, nScenarios=15, eta=.05, liquiditySpread0=.01)



cleanEx()
nameEx("rRealEstate")
### * rRealEstate

flush(stderr()); flush(stdout())

### Name: rRealEstate
### Title: rRealEstate
### Aliases: rRealEstate

### ** Examples

data(ZC)
rRealEstate(horizon=5, nScenarios=10, ZC=ZC, vol=.1, k=2, volRealEstate=.15, realEstate0=50)



cleanEx()
nameEx("rShortRate")
### * rShortRate

flush(stderr()); flush(stdout())

### Name: rShortRate
### Title: rShortRate
### Aliases: rShortRate

### ** Examples

data(ZC)
rShortRate(horizon=15, nScenarios=10, ZC=ZC, vol=.1, k=2)



cleanEx()
nameEx("rStock")
### * rStock

flush(stderr()); flush(stdout())

### Name: rStock
### Title: rStock
### Aliases: rStock

### ** Examples

data(ZC)
rStock(horizon=10, nScenarios=7, ZC=ZC, vol=.1, k=2, volStock=.2, stock0=100, rho=.5)



cleanEx()
nameEx("setForwardRates")
### * setForwardRates

flush(stderr()); flush(stdout())

### Name: setForwardRates
### Title: setForwardRates method
### Aliases: setForwardRates setForwardRates,Scenarios-method

### ** Examples

scenarios1 <- new("Scenarios")
scenarios1 <- setRiskParamsScenarios(scenarios1, vol=.1, k=2,volStock=.2, volRealEstate=.15,volDefault=.2, alpha=.1,beta=1, eta=.05,rho=.5, stock0=100,realEstate0=50, liquiditySpread0=.01, defaultSpread0=.01)
data(ZC)
scenarios1 <- setForwardRates(scenarios1, ZC, horizon=5)



cleanEx()
nameEx("setParamsBaseScenarios")
### * setParamsBaseScenarios

flush(stderr()); flush(stdout())

### Name: setParamsBaseScenarios
### Title: setParamsBaseScenarios method
### Aliases: setParamsBaseScenarios setParamsBaseScenarios,Scenarios-method

### ** Examples

scenarios1 <- new("Scenarios")
scenarios1 <- setParamsBaseScenarios(scenarios1, horizon=5, nScenarios=10)



cleanEx()
nameEx("setRiskParamsScenarios")
### * setRiskParamsScenarios

flush(stderr()); flush(stdout())

### Name: setRiskParamsScenarios
### Title: setRiskParamsScenarios method
### Aliases: setRiskParamsScenarios setRiskParamsScenarios,Scenarios-method

### ** Examples

scenarios1 <- new("Scenarios")
scenarios1 <- setParamsBaseScenarios(scenarios1, horizon=5, nScenarios=10)
scenarios1 <- setRiskParamsScenarios(scenarios1, vol=.1, k=2,volStock=.2, volRealEstate=.15,volDefault=.2, alpha=.1,beta=1, eta=.05,rho=.5, stock0=100,realEstate0=50, liquiditySpread0=.01, defaultSpread0=.01)



cleanEx()
nameEx("setZCRates")
### * setZCRates

flush(stderr()); flush(stdout())

### Name: setZCRates
### Title: setZCRates method
### Aliases: setZCRates setZCRates,Scenarios-method

### ** Examples

scenarios1 <- new("Scenarios")
scenarios1 <- setParamsBaseScenarios(scenarios1, horizon=5, nScenarios=10)
scenarios1 <- setRiskParamsScenarios(scenarios1, vol=.1, k=2,volStock=.2, volRealEstate=.15,volDefault=.2, alpha=.1,beta=1, eta=.05,rho=.5, stock0=100,realEstate0=50, liquiditySpread0=.01, defaultSpread0=.01)
data(ZC)
scenarios1 <- setZCRates(scenarios1, ZC, horizon=5)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
