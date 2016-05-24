########################################################################################################################
## La fonction VarSelModelMLE permet d'effectuer l'estimation des parametres en considerant que les variables donnees
## dans le slot model de l'objet VSLCMresultsContinuous ou VSLCMresultsCategorical.
## Il appelle le code c++ et retourne un objet VSLCMresultsContinuous ou VSLCMresultsCategorical en fonction de la
## nature des donnees.
########################################################################################################################
setGeneric ( name= "VarSelModelMLE",  def = function(obj,it){ standardGeneric("VarSelModelMLE")})
## Pour les variables continues
setMethod( f = "VarSelModelMLE", 
           signature(obj="VSLCMresultsContinuous",it="numeric"), 
           definition = function(obj,it){
             reference <- OptimizeMICL(obj, "Continuous")
             return(reference)         
           }
)

########################################################################################################################
## La fonction VarSelCluster est disponible pour l'utilisateur et permet d'appeller les fonctions VarSelModelSelection et
## VarSelModelMLE et retourne un objet VSLCMresultsContinuous ou VSLCMresultsCategorical
## La fonction possede deux parametres obligatoires:
## x: tableau de donnees sous format data.frame avec variables numeric pour les continues et factor pour les categorielles
## g: nombre de classes (numeric de taille 1)
## La fonction possede egalement 9 parametres optionnels
## vbleSelec: logical indiquant si la selection de variables est effectuee
## paramEstim: logical indiquant si l'estimation des parametres est effectuee
## nbcores: nombre de coeurs de calcul utilises
## nbSmall: nombre d'initialisations du small EM
## iterSmall: nombre d'iterations des small EM
## nbKeep: nombre de chaines conservees apres le small EM
## iterKeep: nombre d'iterations maximum des EM
## tolKeep: difference des vraisemblances de deux iterations successives impliquant un arret de EM
########################################################################################################################
VarSelCluster <- function(x, g, initModel=50, vbleSelec=TRUE, paramEstim=TRUE, nbcores=1, nbSmall=250, iterSmall=20, nbKeep=50, iterKeep=10**3, tolKeep=10**(-3)){
  # Verifie les parametres d'entrees
  discrim <- rep(1,ncol(x))
  CheckInputs(x, g, initModel, vbleSelec, discrim, paramEstim, nbcores, nbSmall, iterSmall, nbKeep, iterKeep, tolKeep)
  # Creation de l'objet S4 VSLCMstrategy contenant les parametres de reglage
  strategy <- VSLCMstrategy(initModel, nbcores, vbleSelec, paramEstim, nbSmall, iterSmall, nbKeep, iterKeep, tolKeep)    
  # Creation de l'objet S4 VSLCMdataContinuous ou VSLCMdataCategorical
  data <- VSLCMdata(x)
  
  if (class(data) == "VSLCMdataContinuous")
    reference <- new("VSLCMresultsContinuous", data=data, criteria=new("VSLCMcriteria", MICL=-Inf), model=new("VSLCMmodel",g=g, omega=discrim), strategy=strategy)
  else
    stop("The data set is not a continuous one")      

  if (g==1){
    reference <- withoutmixture(reference)
  }else{
    # Estimation du modele et/ou des parametres
    if (strategy@parallel == FALSE)
      reference <- VarSelModelMLE(reference, 0)
    else{
      nb.cpus <- min(detectCores(all.tests = FALSE, logical = FALSE) , max(strategy@initModel,1), nbcores)
      if (strategy@vbleSelec == TRUE){
        reference@strategy <- JustModelStrategy(strategy, nb.cpus)
        
        if(Sys.info()["sysname"] == "Windows")
        {
          cl <- makeCluster(nb.cpus)
          common.objects <- c("reference","VarSelModelMLE", "OptimizeMICL")
          clusterEvalQ(cl, {require(VarSelLCM)})
          clusterExport(cl=cl, varlist = common.objects, envir = environment())
          reference <- parLapply(cl = cl, 
                                 X  = as.list(rep(0, nb.cpus)), 
                                 fun = function(g){VarSelModelMLE(reference,g)})
          stopCluster(cl)
          
        }
        else
          reference <- mclapply(X = as.list(rep(0, nb.cpus)),
                                FUN = VarSelModelMLE,
                                obj=reference,
                                mc.cores = nb.cpus, mc.preschedule = TRUE, mc.cleanup = TRUE)
        # On conserve le meilleur modele au sens de MICL
        tmpMICL <- rep(NA, length(reference))
        for (it in 1:length(reference)) tmpMICL[it] <- reference[[it]]@criteria@MICL

        cvrate <- 0
        for (it in which(tmpMICL==max(tmpMICL))){
          cvrate <- cvrate + reference[[it]]@cvrate
        }
          
        
        reference <- reference[[which.max(tmpMICL)]]
        reference@cvrate = cvrate
        # On parallelise aussi pour les EM donc on reparti les initialisations sur les differents coeurs
      }
      reference@strategy <- strategy 
      nb.cpus <- min(detectCores(all.tests = FALSE, logical = FALSE) , max(reference@strategy@nbSmall,1), nbcores, 4)
      if (strategy@paramEstim){
        reference@strategy@vbleSelec <- FALSE
        reference@strategy@nbSmall <- ceiling(reference@strategy@nbSmall / nb.cpus)
        reference@strategy@nbKeep <- ceiling(reference@strategy@nbKeep / nb.cpus)
        
        if(Sys.info()["sysname"] == "Windows")
        {
          cl <- makeCluster(nb.cpus)
          common.objects <- c("reference","VarSelModelMLE", "OptimizeMICL")
          clusterEvalQ(cl, {require(VarSelLCM)})
          clusterExport(cl=cl, varlist = common.objects, envir = environment())
          reference <- parLapply(cl = cl, 
                                 X  = as.list(rep(0, nb.cpus)), 
                                 fun = function(g){VarSelModelMLE(reference,g)})
          stopCluster(cl)
          
        }
        else
          reference <- mclapply(X = as.list(rep(0, nb.cpus)), FUN = VarSelModelMLE, obj=reference, mc.cores = nb.cpus, mc.preschedule = TRUE, mc.cleanup = TRUE)
        
        # On conserve les parametres maximisant la vraisemblance
        tmploglike <- rep(NA, length(reference))
        for (it in 1:length(tmploglike)) {if (reference[[it]]@criteria@degeneracyrate!=1) tmploglike[it] <- reference[[it]]@criteria@loglikelihood}
        if (all(is.na(tmploglike))) tmploglike[1]=1
        reference <- reference[[which.max(tmploglike)]]
        reference@strategy <- strategy
      }
    }
  }
  return(DesignOutput(reference))
}