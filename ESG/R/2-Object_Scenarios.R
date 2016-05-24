#Scenarios object including models' parameters and risk factors paths
#Each path contains nScenarios rows and horizon+1 columns.
#Warning: first data is at time 0, ex: shortRatePaths[,0]=r_0,shortRatePaths[,t+1]=r_t
############################################
#Scenarios class definition



setClass(  
  Class="Scenarios",
  representation=representation(
    ParamsScenarios="ParamsScenarios", #models parameters
    ForwardRates="vector", #initial forward rates for HJM models approach
    ZCRates="vector",
    shortRatePaths="matrix",
    stockPaths="matrix",
    realEstatePaths="matrix",
    liquiditySpreadPaths="matrix",
    defaultSpreadPaths="matrix"
  )
)

############################################
############################################


setGeneric(
  name="setParamsBaseScenarios",
  def = function(.Object, horizon, nScenarios)
  {
    standardGeneric("setParamsBaseScenarios")
  }
)

setMethod(
  f="setParamsBaseScenarios",
  signature="Scenarios",
  definition=function(.Object, horizon, nScenarios)
  {
    .Object@ParamsScenarios@horizon  <- horizon
    .Object@ParamsScenarios@nScenarios   <- nScenarios
    return(.Object)
  }
)


############################################
############################################



setGeneric(
  name="setRiskParamsScenarios",
  def = function(.Object, vol, k,
                 volStock, volRealEstate, volDefault,
                 alpha, beta, eta, rho,
                 stock0, realEstate0, liquiditySpread0,
                 defaultSpread0)
  {
    standardGeneric("setRiskParamsScenarios")
  }
)

setMethod(
  f="setRiskParamsScenarios",
  signature="Scenarios",
  definition=function(.Object, vol, k,
                      volStock, volRealEstate, volDefault,
                      alpha, beta, eta, rho,
                      stock0, realEstate0, liquiditySpread0,
                      defaultSpread0)
  {
    .Object@ParamsScenarios@vol <- vol
    .Object@ParamsScenarios@k	<- k
    .Object@ParamsScenarios@volStock	<- volStock
    .Object@ParamsScenarios@volRealEstate	<- volRealEstate
    .Object@ParamsScenarios@volDefault	<- volDefault
    .Object@ParamsScenarios@alpha	<- alpha
    .Object@ParamsScenarios@beta	<- beta
    .Object@ParamsScenarios@eta	<- eta
    .Object@ParamsScenarios@rho	<- rho
    .Object@ParamsScenarios@stock0	<- stock0
    .Object@ParamsScenarios@realEstate0	<- realEstate0
    .Object@ParamsScenarios@liquiditySpread0	<- liquiditySpread0
    .Object@ParamsScenarios@defaultSpread0	<- defaultSpread0   
    return(.Object)
  }
)

############################################
############################################

setGeneric(
  name="setRiskParamsScenariosrt",
  def = function(.Object, vol, k)
  {
    standardGeneric("setRiskParamsScenariosrt")
  }
)

setMethod(
  f="setRiskParamsScenariosrt",
  signature="Scenarios",
  definition=function(.Object, vol, k)
  {
    .Object@ParamsScenarios@vol <- vol
    .Object@ParamsScenarios@k  <- k
    return(.Object)
  }
)
############################################
############################################


setGeneric(
  name="setRiskParamsScenariosS",
  def = function(.Object, vol, k,
                 volStock, 
                 stock0, rho)
  {
    standardGeneric("setRiskParamsScenariosS")
  }
)

setMethod(
  f="setRiskParamsScenariosS",
  signature="Scenarios",
  definition=function(.Object, vol, k,
                      volStock, 
                      stock0, rho)
  {
    .Object@ParamsScenarios@vol <- vol
    .Object@ParamsScenarios@k  <- k
    .Object@ParamsScenarios@volStock	<- volStock
    .Object@ParamsScenarios@stock0	<- stock0
    .Object@ParamsScenarios@rho  <- rho
    
    return(.Object)
  }
)

############################################
############################################


setGeneric(
  name="setRiskParamsScenariosRE",
  def = function(.Object, vol, k,
                 volRealEstate, 
                 realEstate0)
  {
    standardGeneric("setRiskParamsScenariosRE")
  }
)

setMethod(
  f="setRiskParamsScenariosRE",
  signature="Scenarios",
  definition=function(.Object, vol, k,
                      volRealEstate, 
                      realEstate0)
  {
    .Object@ParamsScenarios@vol <- vol
    .Object@ParamsScenarios@k  <- k
    .Object@ParamsScenarios@volRealEstate  <- volRealEstate
    .Object@ParamsScenarios@realEstate0	<- realEstate0
    return(.Object)
  }
)

############################################
############################################


setGeneric(
  name="setRiskParamsScenariosliqSpr",
  def = function(.Object, eta, liquiditySpread0)
  {
    standardGeneric("setRiskParamsScenariosliqSpr")
  }
)

setMethod(
  f="setRiskParamsScenariosliqSpr",
  signature="Scenarios",
  definition=function(.Object,eta, liquiditySpread0)
  {
    .Object@ParamsScenarios@eta  <- eta
    .Object@ParamsScenarios@liquiditySpread0 <- liquiditySpread0 
    return(.Object)
  }
)

############################################
############################################


setGeneric(
  name="setRiskParamsScenariosdefSpr",
  def = function(.Object, volDefault, defaultSpread0, alpha, beta)
  {
    standardGeneric("setRiskParamsScenariosdefSpr")
  }
)

setMethod(
  f="setRiskParamsScenariosdefSpr",
  signature="Scenarios",
  definition=function(.Object, volDefault, defaultSpread0, alpha, beta)
  {
    .Object@ParamsScenarios@volDefault <- volDefault 
    .Object@ParamsScenarios@defaultSpread0 <- defaultSpread0        
    .Object@ParamsScenarios@alpha  <- alpha
    .Object@ParamsScenarios@beta <- beta
    return(.Object)
  }
)


############################################
############################################



setGeneric(
  name="setForwardRates",
  def = function(.Object, ZC, horizon)
  {
    standardGeneric("setForwardRates")
  }
)

setMethod(
  f="setForwardRates",
  signature="Scenarios",
  definition=function(.Object, ZC, horizon)
  {
    .Object@ForwardRates <- ForwardExtraction(ZC,horizon)    
    return(.Object)
  }
)

############################################
############################################


setGeneric(
  name="setZCRates",
  def = function(.Object, ZC, horizon)
  {
    standardGeneric("setZCRates")
  }
)

setMethod(
  f="setZCRates",
  signature="Scenarios",
  definition=function(.Object, ZC, horizon)
  {
    .Object@ZCRates <- ZCExtraction(ZC, horizon)    
    return(.Object)
  }
)

############################################
############################################



setGeneric(
  name="getParamsBaseScenarios",
  def = function(.Object)
  {
    standardGeneric("getParamsBaseScenarios")
  }
)

setMethod(
  f="getParamsBaseScenarios",
  signature="Scenarios",
  definition=function(.Object)
  {   
    return(list(horizon = .Object@ParamsScenarios@horizon , nScenarios = .Object@ParamsScenarios@nScenarios ))
  }
    )


############################################
############################################



setGeneric(
  name="getRiskParamsScenarios",
  def = function(.Object)

{
  standardGeneric("getRiskParamsScenarios")
}
)

setMethod(
  f="getRiskParamsScenarios",
  signature="Scenarios",
  definition=function(.Object)
  { 
    return(list(vol=.Object@ParamsScenarios@vol, k=.Object@ParamsScenarios@k, 
                volStock=.Object@ParamsScenarios@volStock , volRealEstate =  .Object@ParamsScenarios@volRealEstate,
                volDefault =  .Object@ParamsScenarios@volDefault, alpha= .Object@ParamsScenarios@alpha,
                beta =  .Object@ParamsScenarios@beta,  eta =  .Object@ParamsScenarios@eta, rho =  .Object@ParamsScenarios@rho, stock0 =  .Object@ParamsScenarios@stock0, realEstate0 =  .Object@ParamsScenarios@realEstate0,  liquiditySpread0 =  .Object@ParamsScenarios@liquiditySpread0,  
                defaultSpread0 = .Object@ParamsScenarios@defaultSpread0))
  }
)

############################################
############################################


setGeneric(
  name="getRiskParamsScenariosrt",
  def = function(.Object, vol, k)
  {
    standardGeneric("getRiskParamsScenariosrt")
  }
)

setMethod(
  f="getRiskParamsScenariosrt",
  signature="Scenarios",
  definition=function(.Object)
  {
    return(list(vol = .Object@ParamsScenarios@vol, k = .Object@ParamsScenarios@k))
  }
    )
############################################
############################################



setGeneric(
  name="getRiskParamsScenariosS",
  def = function(.Object)
  {
    standardGeneric("getRiskParamsScenariosS")
  }
)

setMethod(
  f="getRiskParamsScenariosS",
  signature="Scenarios",
  definition=function(.Object)
  {
    return(list(vol=.Object@ParamsScenarios@vol, k=.Object@ParamsScenarios@k, 
                volStock=.Object@ParamsScenarios@volStock, stock0=.Object@ParamsScenarios@stock0))
  }
    )

############################################
############################################


setGeneric(
  name="getRiskParamsScenariosRE",
  def = function(.Object)
  {
    standardGeneric("getRiskParamsScenariosRE")
  }
)

setMethod(
  f="getRiskParamsScenariosRE",
  signature="Scenarios",
  definition=function(.Object)
  {
    return(list(vol=.Object@ParamsScenarios@vol, k=.Object@ParamsScenarios@k, 
                volRealEstate=.Object@ParamsScenarios@volRealEstate, realEstate0 =  .Object@ParamsScenarios@realEstate0))
  }
)

############################################
############################################


setGeneric(
  name="getRiskParamsScenariosliqSpr",
  def = function(.Object)
  {
    standardGeneric("getRiskParamsScenariosliqSpr")
  }
)

setMethod(
  f="getRiskParamsScenariosliqSpr",
  signature="Scenarios",
  definition=function(.Object)
  {
    return(list(eta=.Object@ParamsScenarios@eta, liquiditySpread0=.Object@ParamsScenarios@liquiditySpread0))
  }
)

############################################
############################################


setGeneric(
  name="getRiskParamsScenariosdefSpr",
  def = function(.Object)
  {
    standardGeneric("getRiskParamsScenariosdefSpr")
  }
)

setMethod(
  f="getRiskParamsScenariosdefSpr",
  signature="Scenarios",
  definition=function(.Object)
  {
    return(list(volDefault=.Object@ParamsScenarios@volDefault, alpha=.Object@ParamsScenarios@alpha, 
                beta=.Object@ParamsScenarios@beta,  defaultSpread0=.Object@ParamsScenarios@defaultSpread0))
  }
)


############################################
############################################



setGeneric(
  name="getForwardRates",
  def = function(.Object, ZC, horizon)
  {
    standardGeneric("getForwardRates")
  }
)

setMethod(
  f="getForwardRates",
  signature="Scenarios",
  definition=function(.Object)
  {
    return(.Object@ForwardRates)
  }
)

############################################
############################################


setGeneric(
  name="getZCRates",
  def = function(.Object)
  {
    standardGeneric("getZCRates")
  }
)

setMethod(
  f="getZCRates",
  signature="Scenarios",
  definition=function(.Object)
  {
    return(.Object@ZCRates)
  }
)

############################################
############################################


setGeneric(
  name="getShortRatePaths",
  def = function(.Object)
  {
    standardGeneric("getShortRatePaths")
  }
)

setMethod(
  f="getShortRatePaths",
  signature="Scenarios",
  definition=function(.Object)
  {
    return(.Object@shortRatePaths)
  }
)
  

############################################
############################################


setGeneric(
  name="getstockPaths",
  def = function(.Object)
  {
    standardGeneric("getstockPaths")
  }
)

setMethod(
  f="getstockPaths",
  signature="Scenarios",
  definition=function(.Object)
  {
    return(list(shortRatePaths = .Object@shortRatePaths, stockPaths = .Object@stockPaths))
  }
)
############################################
############################################


setGeneric(
  name="getrealEstatePaths",
  def = function(.Object)
  {
    standardGeneric("getrealEstatePaths")
  }
)

setMethod(
  f="getrealEstatePaths",
  signature="Scenarios",
  definition=function(.Object)
  {
    return(list(shortRatePaths = .Object@shortRatePaths, realEstatePaths = .Object@realEstatePaths))
  }
)
############################################
############################################


setGeneric(
  name="getLiquiditySpreadPaths",
  def = function(.Object)
  {
    standardGeneric("getLiquiditySpreadPaths")
  }
)

setMethod(
  f="getLiquiditySpreadPaths",
  signature="Scenarios",
  definition=function(.Object)
  {
    return(.Object@liquiditySpreadPaths)
  }
)
############################################
############################################


setGeneric(
  name="getdefaultSpreadPaths",
  def = function(.Object)
  {
    standardGeneric("getdefaultSpreadPaths")
  }
)

setMethod(
  f="getdefaultSpreadPaths",
  signature="Scenarios",
  definition=function(.Object)
  {
    return(.Object@defaultSpreadPaths)
  }
)
