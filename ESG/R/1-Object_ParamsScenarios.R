#Classe ParamsScenarios, contenant tous les parametres des modeles
############################################



setClass(
  Class="ParamsScenarios",
  representation=representation(
    horizon = "numeric",
    nScenarios = "numeric",
    #Vasicek Hull&White
    vol="numeric",
    k="numeric",
    #Black & Scholes
    volStock="numeric",
    volRealEstate="numeric",
    stock0="numeric",
    realEstate0="numeric",
    #LMN
    volDefault="numeric",
    alpha="numeric",
    beta="numeric",
    eta="numeric",
    liquiditySpread0="numeric",
    defaultSpread0="numeric",
    #Correlation
    rho="numeric"
  )
)