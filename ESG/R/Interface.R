# User's functions

##### Risk factors 

# Short rate generation


rShortRate <- function(horizon, nScenarios, ZC, vol, k)
{
  # Short rate parameters 
  # ZC : ZC rate input
  # vol, k
  
  # Scenario object creation
  objScenario  <- new("Scenarios") 
  
  # Basic scenario's parameters setting
  objScenario  <- setParamsBaseScenarios(objScenario, horizon = horizon, nScenarios = nScenarios)
  
  # Risk factors parameters setting
  objScenario  <- setRiskParamsScenariosrt(objScenario, vol = vol, k = k)
  
  # Forward and ZC rates setting
  objScenario  <- setForwardRates(objScenario, ZC, horizon) 
  objScenario  <- setZCRates(objScenario, ZC, horizon)
  
  # Projection
  objScenario  <- customPathsGeneration(objScenario, type="shortRate")
  
  y <- getShortRatePaths(objScenario)   
  
  # Object deletation
  rm(objScenario)
  
  return (y)
}

# Stock generation



rStock <- function(horizon, nScenarios, ZC, vol, k, volStock, stock0, rho)
{
  # Short rate parameters 
  # ZC : ZC rate input
  # vol, k
  # Stock parameters
  # volStock, stock0
  
  # Scenario object creation
  objScenario  <- new("Scenarios") 
  
  # Basic scenario's parameters setting
  objScenario  <- setParamsBaseScenarios(objScenario, horizon = horizon, nScenarios = nScenarios)
  
  # Risk factors parameters setting
  objScenario  <- setRiskParamsScenariosrt(objScenario, vol = vol, k = k)
  objScenario  <- setRiskParamsScenariosS(objScenario, vol = vol, k = k, volStock = volStock, stock0 = stock0, rho=rho)
  
  # Forward and ZC rates setting
  objScenario  <- setForwardRates(objScenario, ZC, horizon) 
  objScenario  <- setZCRates(objScenario, ZC, horizon)
  
  # Projection
  objScenario  <- customPathsGeneration(objScenario, type="shortRate")
  objScenario  <- customPathsGeneration(objScenario, type="stock")
  y <- getstockPaths(objScenario)   
  
  # Scenario deletion
  rm(objScenario)
  
  return (y)
}


# Real-Estate generation
rRealEstate <- function(horizon, nScenarios, ZC, vol, k, volRealEstate, realEstate0)
{
  # Short rate parameters 
  # ZC : ZC rate input
  # vol, k
  # Real-estate parameters
  # volRealEstate, realEstate0
  
  # Scenario object creation
  objScenario  <- new("Scenarios") 
  
  # Basic scenario's parameters setting
  objScenario  <- setParamsBaseScenarios(objScenario, horizon = horizon, nScenarios = nScenarios)
  
  # Risk factors parameters setting
  objScenario  <- setRiskParamsScenariosrt(objScenario, vol = vol, k = k)
  objScenario  <- setRiskParamsScenariosRE(objScenario, vol = vol, k = k, volRealEstate = volRealEstate, realEstate0 = realEstate0)
  
  # Forward and ZC rates setting
  objScenario  <- setForwardRates(objScenario, ZC, horizon) 
  objScenario  <- setZCRates(objScenario, ZC, horizon)

  # Projection
  objScenario  <- customPathsGeneration(objScenario, type="shortRate")
  objScenario  <- customPathsGeneration(objScenario, type="realEstate")
  
  y <- getrealEstatePaths(objScenario)
  
  # Deletion
  rm(objScenario)
  
  return (y)
}




# Liquidity spread generation
rLiquiditySpread <- function(horizon, nScenarios, eta, liquiditySpread0)
{
  # Liquidity spread parameters 
  # eta, liquiditySpread0
  
  # Scenario object creation
  objScenario  <- new("Scenarios") 
  
  # Basic parameters setting
  objScenario  <- setParamsBaseScenarios(objScenario, horizon = horizon, nScenarios = nScenarios)
  
  # Risk factor parameters setting
  objScenario  <- setRiskParamsScenariosliqSpr(objScenario, eta, liquiditySpread0)
  
  objScenario  <- customPathsGeneration(objScenario, type="liquiditySpread")
  
  y <- getLiquiditySpreadPaths(objScenario)
  
  # Object deletion
  rm(objScenario)
  
  return(y)
}


# Default spread generation


rDefaultSpread <- function(horizon, nScenarios, defaultSpread0, volDefault, alpha, beta)
{
  # Spread generation 
  # defaultSpread0, volDefault, alpha, beta
  
  # Scenario object creation
  objScenario  <- new("Scenarios") 
  
  # Basic parameters setting
  objScenario  <- setParamsBaseScenarios(objScenario, horizon = horizon, nScenarios = nScenarios)
  
  # Risk factor parameters setting
  objScenario  <-setRiskParamsScenariosdefSpr(objScenario, volDefault, defaultSpread0, alpha, beta)
  
  # Projection
  objScenario  <- customPathsGeneration(objScenario, type="defaultSpread")
  
  y <- getdefaultSpreadPaths(objScenario)
  
  # Object deletion
  rm(objScenario)
  
  return(y)
}




# All risk factors generation 
rAllRisksFactors <- function(horizon, nScenarios, ZC, vol, k, volStock, stock0, rho, volRealEstate, realEstate0, 
                             eta, liquiditySpread0, defaultSpread0, volDefault, alpha, beta)
{ 
  # Short rate parameters 
  # ZC : Zero coupon rates
  # vol, k
  # Stock parameters
  # volStock, stock0
  # Real-Estate parameters
  # volRealEstate, realEstate0
  # Liquidity spread parameters 
  # eta, liquiditySpread0  
  # Default spread parameters
  # defaultSpread0, volDefault, alpha, beta
  
  # Scenario object creation
  objScenario  <- new("Scenarios") 
  
  # Basic parameters setting
  objScenario  <- setParamsBaseScenarios(objScenario, horizon = horizon, nScenarios = nScenarios)
  
  # Forward and zero coupon rates setting
  objScenario  <- setForwardRates(objScenario, ZC, horizon) 
  objScenario  <- setZCRates(objScenario, ZC, horizon)
  
  # Risk factors parameters setting
  # Short rate and stock
  objScenario  <- setRiskParamsScenariosS(objScenario, vol = vol, k = k, volStock = volStock, stock0 = stock0, rho=rho)
  # Real-Estate
  objScenario  <- setRiskParamsScenariosRE(objScenario, vol = vol, k = k, volRealEstate = volRealEstate, realEstate0 = realEstate0)
  # Liquidity
  objScenario  <- setRiskParamsScenariosliqSpr(objScenario, eta=eta, liquiditySpread0=liquiditySpread0)
  # Default
  objScenario  <- setRiskParamsScenariosdefSpr(objScenario, volDefault=volDefault, defaultSpread0=defaultSpread0, alpha=alpha, beta=beta)
  
  # Path generation
  objScenario  <- customPathsGeneration(objScenario)
  
  y <- list(shortRate=getShortRatePaths(objScenario),
            s=getstockPaths(objScenario)$stock,
            realEstate=getrealEstatePaths(objScenario)$realEstate,
            liquiditySpread=getLiquiditySpreadPaths(objScenario),
            defaultSpread=getdefaultSpreadPaths(objScenario))
  
  # Object deletion
  rm(objScenario)
  
  return(y)
}

##### Assets values


rAssetDistribution <- function(type, t, T, vol, k, ZC, nScenarios = NULL, 
                               volStock=NULL, stock0=NULL, rho=NULL, volRealEstate=NULL, realEstate0=NULL, 
                               eta=NULL, liquiditySpread0=NULL, defaultSpread0=NULL, volDefault=NULL, alpha=NULL, beta=NULL,
                               nCoupons=NULL,couponsRate=NULL,omega=NULL, s=NULL, Strike=NULL)
{ 
  if (missing(type)||!(type %in% c("Zero-Coupon", "Bond", "CBond", "ConvBond", "EuroCall_S", "EuroPut_Stock", "EuroCall_ZC", "EuroPut_ZC", "CDS"))) 
  {stop("Parameter 'type' must be : 'Zero-Coupon', 'Bond', 'CBond', 'ConvBond', 'EuroCall_S', 'EuroPut_Stock', 'EuroCall_ZC', 'EuroPut_ZC', 'CDS' ")}
  
  if (missing(t)||missing(T)) {stop("You must fill at least the parameters t and T.")} else
  {
    ########## Verifs    
    if (t == 0 && missing(nScenarios)) {nScenarios <- 1}
    
    if (t == 0 && nScenarios>1) 
    {
      nScenarios <- 1
      warning("When t = 0, parameter nScenarios is unused and set to 1")
    }
    
    if (t > 0 && (missing(nScenarios)||nScenarios<=1)) {stop("When parameter t > 0, parameter nScenarios should be provided and >= 2")}
    
    if (t > 0 && missing(ZC)){stop("Zero-coupon yield curve (parameter ZC) is missing")} 
    
    if (t > 0 && (!missing(ZC)))
    {
      maxT  <- length(ZC)/12
      if(t > maxT||T > maxT){stop("Not enough yield curve maturities found for this calculation. For the given yield curve, both starting month t=",t," and maturity T=",T," must be <=",floor(maxT))}  
    }       
    ##########
    
    # Scenario object creation
    objScenario  <- new("Scenarios") 
    
    # Basic parameters setting
    horizon <- ceiling(max(T,s))
    objScenario  <- setParamsBaseScenarios(objScenario, horizon = horizon, nScenarios = nScenarios)
    
    # Forward and zero-coupon rates
    objScenario  <- setForwardRates(objScenario, ZC, horizon) 
    objScenario  <- setZCRates(objScenario, ZC, horizon)
    objScenario  <- setRiskParamsScenariosrt(objScenario, vol = vol, k = k)      
    
    if (type == "Zero-Coupon")  
    {
      if (!missing(nCoupons)||!missing(couponsRate)||!missing(omega)||!missing(s)||!missing(Strike)) {warning("Only these parameters are actually used for this asset : t, T")}  
      objScenario  <- customPathsGeneration(objScenario, type="shortRate")
      y <- Asset_PriceDistribution(objScenario,type=type,t=t,T=T)
    }
    
    if (type == "Bond") 
    {    
      if (missing(nCoupons)||missing(couponsRate)) {stop("You must fill the parameters nCoupons and couponsRate.")} else {
        if (!missing(omega)||!missing(s)||!missing(Strike)) {warning("Only these parameters are actually used for this asset : t, T, nCoupons, couponsRate")}
        objScenario  <- customPathsGeneration(objScenario, type="shortRate")
        y <- Asset_PriceDistribution(objScenario,type=type,t=t,T=T,nCoupons=nCoupons,couponsRate=couponsRate)}
    }
    
    if (type == "CBond") 
    {
      if (missing(nCoupons)||missing(couponsRate)||missing(omega)) {stop("You must fill the parameters nCoupons, couponsRate and omega.")} else {
        if (!missing(s)||!missing(Strike)) {warning("Only these parameters are actually used for this asset : t,T,nCoupons,couponsRate,omega")}
        if(missing(eta)||missing(liquiditySpread0)||missing(volDefault)||missing(defaultSpread0)||missing(alpha)||missing(beta))
        {stop("Parameters eta, liquiditySpread0, volDefault, defaultSpread0, alpha and beta should be provided")}
        
        objScenario  <- customPathsGeneration(objScenario, type="shortRate")
        objScenario  <- setRiskParamsScenariosliqSpr(objScenario, eta=eta, liquiditySpread0=liquiditySpread0)
        objScenario  <- setRiskParamsScenariosdefSpr(objScenario, volDefault=volDefault, defaultSpread0=defaultSpread0, alpha=alpha, beta=beta)        
        objScenario  <- customPathsGeneration(objScenario, type="liquiditySpread")
        objScenario  <- customPathsGeneration(objScenario, type="defaultSpread")
        y <- Asset_PriceDistribution(objScenario,type=type,t=t,T=T,nCoupons=nCoupons,couponsRate=couponsRate,omega=omega)}
      
    }
    
    if (type == "EuroCall_Stock")  
    {
      
      if (missing(Strike)) {stop("You must fill the Strike parameter")} else {
        if (!missing(nCoupons)||!missing(couponsRate)||!missing(omega)||!missing(s)) {warning("Only these parameters are actually used for this asset: t, T, Strike")}
        if (missing(volStock)||missing(stock0)||missing(rho)){stop("Parameters volStock, stock0, rho should be provided")}
        objScenario  <- setRiskParamsScenariosS(objScenario, vol = vol, k = k, volStock = volStock, stock0 = stock0, rho=rho)        
        objScenario  <- customPathsGeneration(objScenario, type="stock")
        y <- Asset_PriceDistribution(objScenario,type=type,t=t,T=T,Strike=Strike)}
    }
    
    if (type == "EuroPut_Stock")  
    {
      if (missing(Strike)) {stop("You must fill the Strike parameter")} else {  
        if (!missing(nCoupons)||!missing(couponsRate)||!missing(omega)||!missing(s)) {warning("Only these parameters are actually used for this asset : t, T, Strike")}
        if (missing(volStock)||missing(stock0)||missing(rho)){stop("Parameters volStock, stock0, rho should be provided")}
        objScenario  <- setRiskParamsScenariosS(objScenario, vol = vol, k = k, volStock = volStock, stock0 = stock0, rho=rho)        
        objScenario  <- customPathsGeneration(objScenario, type="stock")
        y <- Asset_PriceDistribution(objScenario,type=type,t=t,T=T,Strike=Strike)}
    }
    
    if (type == "EuroCall_ZC")  
    { 
      if (missing(s)||missing(Strike)) {stop("You must fill the s and Strike parameters")} else {
        if (!missing(nCoupons)||!missing(couponsRate)||!missing(omega)) {warning("Only these parameters are actually used for this asset : t, T, s, Strike")}
        objScenario  <- setRiskParamsScenariosrt(objScenario, vol = vol, k = k)
        objScenario  <- customPathsGeneration(objScenario, type="shortRate")
        y <- Asset_PriceDistribution(objScenario,type=type,t=t,T=T,s=s,Strike=Strike)}
    }
    
    if (type == "EuroPut_ZC")  
    {
      if (missing(s)||missing(Strike)) {stop("You must fill the s and Strike parameters")} else {
        if (!missing(nCoupons)||!missing(couponsRate)||!missing(omega)) {warning("Only these parameters are actually used for this asset : t, T, s, Strike")}
        objScenario  <- setRiskParamsScenariosrt(objScenario, vol = vol, k = k)
        objScenario  <- customPathsGeneration(objScenario, type="shortRate")
        y <- Asset_PriceDistribution(objScenario,type=type,t=t,T=T,s=s,Strike=Strike)}
    }
    
    if (type == "CDS")
    {
      if (missing(omega)||missing(t)||missing(T))
      {
        stop("You must fill the omega, t and T parameters")
      }
      else
      {
        objScenario  <- setRiskParamsScenariosliqSpr(objScenario, eta=eta, liquiditySpread0=liquiditySpread0)
        objScenario  <- setRiskParamsScenariosdefSpr(objScenario, volDefault=volDefault, defaultSpread0=defaultSpread0, alpha=alpha, beta=beta)        
        objScenario  <- customPathsGeneration(objScenario, type="shortRate")
        objScenario  <- customPathsGeneration(objScenario, type="liquiditySpread")
        objScenario  <- customPathsGeneration(objScenario, type="defaultSpread")
        y <- Asset_PriceDistribution(objScenario,type=type,omega=omega,t=t,T=T)}
    }
    
    if (type == "ConvBond")
    {
      if (missing(nCoupons)||missing(couponsRate))
      {
        stop("You must fill the omega, t and T parameters")
      }
      else
      {
        objScenario  <- setRiskParamsScenariosS(objScenario, vol = vol, k = k, volStock = volStock, stock0 = stock0, rho=rho)        
        objScenario  <- customPathsGeneration(objScenario, type="shortRate")
        objScenario  <- customPathsGeneration(objScenario, type="stock")
        y <- Asset_PriceDistribution(objScenario,type=type,t=t,T=T,nCoupons=nCoupons,couponsRate=couponsRate)}
    }
  }
  
  rm(objScenario)
  
  return (y)
  
}