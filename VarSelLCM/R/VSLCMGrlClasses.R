########################################################################################################################
## Differentes classes S4 accessibles a l'utilisateur et etant des slots des classes S4
## VSLCMresultsContinuous et/ou VSLCMresultsCategorical
########################################################################################################################

########################################################################################################################
## Classe S4 VSLCMcriteria contenant la logvraisemblance (loglikelihood), la valeur des criteres BIC, ICL et MICL
########################################################################################################################
setClass(Class = "VSLCMcriteria", 
         representation = representation(loglikelihood="numeric", BIC="numeric", ICL="numeric", MICL="numeric", degeneracyrate="numeric"), 
         prototype = prototype(loglikelihood=numeric(), BIC=numeric(), ICL=numeric(), MICL=numeric(), degeneracyrate=numeric() )
)

########################################################################################################################
## Classe S4 VSLCMpartitions contenant la partition MAP (zMAP), la partition zstar (zOPT) et la partition floue (tik)
########################################################################################################################
setClass(
  Class = "VSLCMpartitions", 
  representation = representation(zMAP="numeric", zOPT="numeric", tik="matrix"), 
  prototype = prototype(zMAP=numeric(), zOPT=numeric(), tik=matrix(0,0,0))
)

########################################################################################################################
## Classe S4 VSLCMstrategy contenant les parametres de reglages detailles dans VarSELLCMmixte.R
########################################################################################################################
setClass(
  Class = "VSLCMstrategy", 
  representation = representation(initModel="numeric", vbleSelec="logical", paramEstim="logical", parallel="logical",
    nbSmall="numeric", iterSmall="numeric", nbKeep="numeric", iterKeep="numeric", tolKeep="numeric"), 
  prototype = prototype(initModel=numeric(), vbleSelec=logical(), paramEstim=logical(), parallel=logical(),
    nbSmall=numeric(), iterSmall=numeric(), nbKeep=numeric(), iterKeep=numeric(), tolKeep=numeric())
) 

## Constructeur de la classe S4 VSLCMstrategy
VSLCMstrategy <- function(initModel=50, nbcores=1, vbleSelec=TRUE, paramEstim=TRUE, nbSmall=100, iterSmall=20, nbKeep=10, iterKeep=100, tolKeep=0.001){
  if( nbKeep > nbSmall)
    nbKeep <- nbSmall
  new("VSLCMstrategy",initModel=initModel, parallel=(nbcores>1), vbleSelec=vbleSelec, paramEstim=paramEstim, nbSmall=nbSmall, iterSmall=iterSmall, nbKeep=nbKeep,
      iterKeep=iterKeep, tolKeep=tolKeep)
}

## Utilise lors de la parallellisation
JustModelStrategy <- function(strategy, nb.cpus){
  output <- strategy
  output@paramEstim <- FALSE
  output@initModel <- ceiling(strategy@initModel / nb.cpus)
  return(output)
}

########################################################################################################################
## Classe S4 VSLCMmodel contenant le nombre de classes (g) et le role des variables (omega)
########################################################################################################################
setClass(
  Class = "VSLCMmodel", 
  representation = representation(g="numeric", omega="numeric"), 
  prototype = prototype(g=numeric(), omega=numeric())
)

########################################################################################################################
## Classe S4 VSLCMparamContinuous contenant les proportions (pi), les moyennes (mu) et les ecrat-types (sd)
########################################################################################################################
setClass(
  Class = "VSLCMparamContinuous", 
  representation = representation(pi="numeric", mu="matrix", sd="matrix"), 
  prototype = prototype(pi=numeric(), mu=matrix(), sd=matrix())
)


########################################################################################################################
## Classe S4 VSLCMresultsContinuous
########################################################################################################################
setClass(
  Class = "VSLCMresultsContinuous", 
  representation = representation(data="VSLCMdataContinuous", criteria="VSLCMcriteria", partitions="VSLCMpartitions",
                                  model="VSLCMmodel", strategy="VSLCMstrategy", param="VSLCMparamContinuous", cvrate="numeric"), 
  prototype = prototype(data=new("VSLCMdataContinuous"), criteria=new("VSLCMcriteria"), partitions=new("VSLCMpartitions"),
                        model=new("VSLCMmodel"), strategy=new("VSLCMstrategy"), param=new("VSLCMparamContinuous"), cvrate=numeric())
)