##' @name BIOMOD_tuning
##' @aliases BIOMOD_tuning
##' 
##' @title Tune models parameters
##' @description Function to tune biomod single models parameters
##'
##' @param data            BIOMOD.formated.data object returned by BIOMOD_FormatingData
##' @param models          vector of models names choosen among 'GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'FDA', 'MARS', 'RF' and 'MAXENT.Phillips'
##' @param models.options  BIOMOD.models.options object returned by BIOMOD_ModelingOptions. Default: BIOMOD_ModelingOptions()
##' @param method.RF       which classification or regression model to use for randomForest (default: "rf"). see http://topepo.github.io/caret/Random_Forest.html
##' @param method.ANN      which classification or regression model to use for artificial neural networks (default: "avNNet"). see http://topepo.github.io/caret/Neural_Network.html
##' @param method.MARS     which classification or regression model to use for mars (default: "earth"). see http://topepo.github.io/caret/Multivariate_Adaptive_Regression_Splines.html
##' @param method.GAM      which classification or regression model to use for GAM (default: "gam"). see http://topepo.github.io/caret/Generalized_Additive_Model.html
##' @param method.GLM      which classification or regression model to use for GLM: (default: 'glmStepAIC'). see http://topepo.github.io/caret/Generalized_Linear_Model.html
##' @param type.GLM        either 'simple', 'quadratic', 'polynomial' or 's_smoother' defining the type of formula you want to build with GLM
##' @param cvmethod.ME     method used for data partitioning for MAXENT.Phillips (default: 'randomkfold')
##' @param kfolds.ME       number of bins to use for k-fold cross-validation used for MAXENT.Phillips (Default: 10).
##' @param overlap.ME      logical; Calculates pairwise metric of niche overlap if TRUE (Default: FALSE). (see ?calc.niche.overlap)
##' @param clamp.ME        logical; If TRUE (Default) "Features are constrained to remain within the range of values in the training data" (Elith et al. 2011)
##' @param n.bg.ME         Number of Background points used to run MAXENT.Phillips (Default: 10000)
##' @param env.ME          RasterStack of model predictor variables
##' @param bin.output.ME   logical; If TRUE, appends evaluations metrics for each evaluation bin to results table (i.e., in addition to the average values across bins).
##' @param trControl       global control parameters for runing (default trainControl(method="cv",summaryFunction=twoClassSummary,classProbs=T),returnData = FALSE). for details see trainControl
##' @param ctrl.GAM        specify control parameters only for GAM (default trControl)
##' @param ctrl.GBM        specify control parameters only for GBM (default trControl)
##' @param ctrl.GLM        specify control parameters only for GLM (default trControl)
##' @param ctrl.CTA        specify control parameters only for CTA (default trControl)
##' @param ctrl.RF         specify control parameters only for RF (default trControl)
##' @param ctrl.ANN        specify control parameters only for ANN (default trControl)
##' @param ctrl.MARS       specify control parameters only for MARS (default trControl)
##' @param ctrl.FDA        specify control parameters only for FDA (default trControl)
##' @param metric          metric to select the optimal model (Default ROC). TSS (maximizing Sensitivity and Specificity) is also possible. see ?train
##' @param metric.ME       metric to select the optimal model for Maxent (Default: ROC). One out of Mean.AUC (or ROC), delta.AICc, Mean.AUC.DIFF. see ?ENMevaluate 
##' @param tuneLength      see ?train (default 30)
##' @param size.tune.ANN   size parameters (number of units in the hidden layer) for ANN used for tuning (default: c(2,4,6,8)).  Will be optimised using the method specified in ctrl.ANN (if not available trControl).
##' @param decay.tune.ANN  weight decay parameters used for tuning for ANN (default: c(0.001, 0.01, 0.05, 0.1))  Will be optimised by method specified in ctrl.ANN (if not available trControl).
##' @param maxit.ANN       maximum number of iterations for ANN (default 500) 
##' @param MaxNWts.ANN     The maximum allowable number of weights for ANN (default 10 * (ncol(myBiomodData'at'data.env.var) + 1) + 10 + 1). 
##  at param parallel        logical. If TRUE, the parallel computing is enabled (highly recommended)
##  at param number.cores    numeric. Number of cores used for parallel computing. (Default: detectCores()-1)
##' 
##' @details
##' Tuning ANN: if no parameters are specified for ANN, they are tuned internally
##'   in BIOMOD_Modelling using cross-validation. ANN parameters are therefore 
##'   tuned in every case either within ecospat.BIOMOD.tunning or BIOMOD_Modelling 
##' 
##' @return
##' BIOMOD.models.options object with optimized parameters
##' 
##' @author Frank Breiner
##' 
##' @seealso BIOMOD_ModelingOptions(), train
##' 
##' @examples
##' \dontrun{
##' # species occurrences
##' DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
##'                                     package="biomod2"))
##' head(DataSpecies)
##' 
##' # the name of studied species
##' myRespName <- 'GuloGulo'
##' 
##' # the presence/absences data for our species 
##' myResp <- as.numeric(DataSpecies[,myRespName])
##' 
##' # the XY coordinates of species data
##' myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]
##' 
##' # Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' myExpl = stack( system.file( "external/bioclim/current/bio3.grd", 
##'                              package="biomod2"),
##'                 system.file( "external/bioclim/current/bio4.grd", 
##'                              package="biomod2"), 
##'                 system.file( "external/bioclim/current/bio7.grd", 
##'                              package="biomod2"),  
##'                 system.file( "external/bioclim/current/bio11.grd", 
##'                              package="biomod2"), 
##'                 system.file( "external/bioclim/current/bio12.grd", 
##'                              package="biomod2"))
##' # 1. Formatting Data
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' 
##' # 2. Defining Models Options using default options.
##' ### Duration for turing all models sequential with default settings 
##' ### on 3.4 GHz processor: approx. 45 min tuning all models in parallel
##' ### (on 8 cores) using foreach loops runs much faster: approx. 14 min
##' 
##' #library(doParallel);cl<-makeCluster(8);registerDoParallel(cl) 
##' 
##' 
##' time.seq<-system.time(Biomod.tuning <- BIOMOD_tuning(myBiomodData,
##'                                                              env.ME = myExpl,
##'                                                              n.bg.ME = ncell(myExpl)))
##' #stopCluster(cl)
##' 
##' myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
##'                                      models = c('RF','CTA'), 
##'                                      models.options = Biomod.tuning$models.options, 
##'                                      NbRunEval=1, 
##'                                      DataSplit=100, 
##'                                      VarImport=0, 
##'                                      models.eval.meth = c('ROC'),
##'                                      do.full.models=FALSE,
##'                                      modeling.id="test")
##' 
##' 
##' #  eval.plot(Biomod.tuning$tune.MAXENT.Phillips at results)
##' par(mfrow=c(1,3))
##' plot(Biomod.tuning$tune.CTA.rpart)
##' plot(Biomod.tuning$tune.CTA.rpart2)
##' plot(Biomod.tuning$tune.RF)
##' }
BIOMOD_tuning <- function(data,
                                  models = c('GLM','GBM','GAM','CTA','ANN','FDA','MARS','RF','MAXENT.Phillips'),
                                  models.options = BIOMOD_ModelingOptions(),
                                  method.ANN = 'avNNet',
                                  method.RF = 'rf',
                                  method.MARS = 'earth',
                                  method.GAM = 'gam',
                                  method.GLM = 'glmStepAIC',
                                  trControl = NULL,
                                  metric = 'ROC',
                                  ctrl.CTA = NULL,
                                  ctrl.RF = NULL,
                                  ctrl.ANN = NULL,
                                  ctrl.MARS = NULL,                                  
                                  ctrl.FDA = NULL,
                                  ctrl.GAM = NULL,
                                  ctrl.GBM = NULL,
                                  ctrl.GLM = NULL,
                                  tuneLength = 30,
                                  decay.tune.ANN = c(0.001, 0.01, 0.05, 0.1),
                                  size.tune.ANN = c(2,4,6,8),
                                  maxit.ANN = 500,
                                  MaxNWts.ANN = 10 * (ncol(data@data.env.var) + 1) + 10 + 1,
                                  type.GLM = 'simple',
                                  cvmethod.ME = 'randomkfold',
                                  overlap.ME = FALSE,
                                  bin.output.ME =TRUE,
                                  kfolds.ME = 10,
                                  n.bg.ME = 10000,
                                  env.ME = NULL,
                                  metric.ME = 'ROC',
                                  clamp.ME = T){

  ## MAXENT.Phillips: http://cran.r-project.org/web/packages/ENMeval/ENMeval.pdf --> ENMevaluate()
  ## or:    http://cran.r-project.org/web/packages/maxent/maxent.pdf -->  tune.maxent()
  packages <- NULL
  if(sum(c('GLM','GBM','GAM','CTA','ANN','FDA','MARS','RF','MAXENT.Phillips') %in% models)>0){
    if(!isNamespaceLoaded("caret")){requireNamespace("caret", quietly = TRUE)}
    if(is.null(trControl)){
      trControl <- caret::trainControl(method = "cv", 
                                        summaryFunction = caret::twoClassSummary,
                                        classProbs = T, 
                                        returnData = F)
    }
  }
  
  if("MAXENT.Phillips" %in% models){
    if(is.null(env.ME)){stop("env.ME argument required!")
    }
  }  
  
  tune.GLM <- tune.MAXENT.Phillips <- tune.GAM <- tune.GBM <- tune.CTA.rpart <- tune.CTA.rpart2 <- tune.RF <- tune.ANN <- tune.MARS <- tune.FDA <- NULL
if('SRE' %in% models){cat("No tuning for SRE!")} 

resp <- data@data.species
if(metric == 'ROC' | metric == 'TSS'){resp <- as.factor(ifelse(resp == 1, "Presence", "Absence"))}

if('GBM' %in% models){  
  
  if(is.null(ctrl.GBM)){ctrl.GBM <- trControl}
  cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n",
            "Start tuning GBM. Start coarse tuning\n"))
  
  tune.grid <- expand.grid(.interaction.depth = seq(2, 8, by = 3),
                           .n.trees = c(500, 1000, 2500),
                           .shrinkage = c(0.001, 0.01, 0.1),
                           .n.minobsinnode = 10)
  
  try(tune.GBM <- caret::train(data@data.env.var, resp,
                        method = "gbm",
                        tuneGrid = tune.grid,
                        trControl = ctrl.GBM,
                        verbose = FALSE))
  cat("Best optimization of coarse tuning:\n")
  cat(paste(tune.GBM$bestTune,"\n-=-=-=-=-=-=-=-=-=-=\n"))
  
  if(!is.null(tune.GBM)){
    cat("Start fine tuning\n")
    
    if(tune.GBM$bestTune$n.trees==2500){
      cat("Best optimization with large trees! Tuning GBM will take a while (approx > 30min).\n")
      n.trees <- seq(2500, 10000, by = 2500)}
    
    if(tune.GBM$bestTune$n.trees==1000){n.trees <- seq(750, 3000, by = 250)}
    
    if(tune.GBM$bestTune$n.trees==500){n.trees <- seq(100, 1000, by = 50)}
    
    tune.grid <- expand.grid(.interaction.depth = c(tune.GBM$bestTune$interaction.depth-1,tune.GBM$bestTune$interaction.depth,tune.GBM$bestTune$interaction.depth+1),
                             .n.trees = n.trees,
                             .shrinkage = c(tune.GBM$bestTune$shrinkage/2,tune.GBM$bestTune$shrinkage,tune.GBM$bestTune$shrinkage*5),
                             .n.minobsinnode = 10)
    tune.GBM <- NULL
    try(tune.GBM <- caret::train(data@data.env.var, resp,
                          method = "gbm",
                          tuneGrid = tune.grid,
                          trControl = ctrl.GBM,
                          verbose = FALSE))  
    }
  cat(paste("\n Finished tuning GBM\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  if('RF' %in% models){
    cat("Start tuning RF\n")
    
    if(is.null(ctrl.RF)){ctrl.RF <- trControl}
    ## give both mtry as bestTune
    try(tune.RF <- caret::train(data@data.env.var, resp,
                         method = method.RF,
                         tuneLength = tuneLength,
                         trControl = ctrl.RF,
                         metric = metric))
    cat(paste("Finished tuning RF\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  if('ANN' %in% models){
    cat("Start tuning ANN\n")
    
    if(is.null(ctrl.ANN)){ctrl.ANN <- trControl}
    ## already tuning: 
    # size: optimised by cross validation based on model AUC (NbCv cross validation; tested size will be the following c(2,4,6, 8))
    # decay: optimised by cross validation on model AUC (NbCv cross validation; tested decay will be the following c(0.001, 0.01, 0.05, 0.1)).
    # could increase maxit from 200 to 500
    # a nice option would be to use model averaging for ann: avNNet in package(caret)
    
    ## Box-cox transformation of predictors
    #box.cox.trans <- data@data.env.var
    ## if Some have zero values, we need to add one to them so that
    ## we can use the Box-Cox transformation. Alternatively, we could
    ## use the Yeo-Johnson transformation without altering the data.   
    #if(sum(box.cox.trans)>0){box.cox.trans <- box.cox.trans + 1}
    #pp <- preProcess(data@data.env.var, method = "BoxCox")
    #contPredTrain <- predict(pp, box.cox.trans)
    
    ## Create a specific candidate set of models to evaluate:
    tune.grid <- expand.grid(.decay = decay.tune.ANN,
                             .size = size.tune.ANN,
                             ## The next option is to use bagging (see the
                             ## next chapter) instead of different random
                             ## seeds.
                             .bag = FALSE)
    try(tune.ANN <- caret::train(data@data.env.var, resp, 
                          method = method.ANN,
                          tuneGrid = tune.grid,
                          trControl = ctrl.ANN,
                          ## Automatically standardize data prior to modeling
                          ## and prediction
                          preProc = c("center", "scale"),
                          linout = TRUE,
                          trace = FALSE,
                          MaxNWts.ANN = MaxNWts.ANN,
                          maxit = maxit.ANN,
                          metric = metric))
    cat(paste("Finished tuning ANN\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  if('GAM' %in% models){
    cat("Start tuning GAM\n")
    
    if(is.null(ctrl.GAM)){ctrl.GAM <- trControl}
    
    try(tune.GAM <-   caret::train(data@data.env.var, resp, 
                            method = method.GAM,
                            trControl = ctrl.GAM))
    cat(paste("Finished tuning GAM\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  
  if('MAXENT.Phillips' %in% models){
    cat("Start tuning MAXENT.Phillips\n")
    cat("Switch down so far..")
#     if(cvmethod.ME != 'randomkfold'){kfolds.ME <- NA}
#     try(tune.MAXENT.Phillips <- ENMevaluate(data@coord[data@data.species==1,],env.ME,bg.coords= data@coord[data@data.species==0,],
#                                    method=cvmethod.ME, kfolds = kfolds.ME, overlap=overlap.ME, 
#                                    bin.output=TRUE, clamp=clamp.ME))
    cat(paste("Finished tuning MAXENT.Phillips\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  if('MARS' %in% models){
    cat("Start tuning MARS\n")
    
    if(is.null(ctrl.MARS)){ctrl.MARS <- trControl}
    
    tune.grid <- expand.grid(.degree = 1:2, .nprune = 2:38)
    try(tune.MARS <-   caret::train(data@data.env.var, resp, 
                             method = method.MARS,
                             tuneGrid = tune.grid,
                             trControl = ctrl.MARS))
    cat(paste("Finished tuning MARS\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  
  if('GLM' %in% models){
    cat("Start tuning GLM\n")
    
    if(is.null(ctrl.GLM)){ctrl.GLM <- trControl}
    try(tune.GLM <-   caret::train( makeFormula("resp",data@data.env.var, type= type.GLM,interaction.level = 0),
                             data=cbind(data@data.env.var,resp=resp),
                             method = method.GLM,
                             trControl = ctrl.GLM))
    cat(paste("Finished tuning GLM\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  
  
  if('FDA' %in% models){
    ## DOES NOT WORK
    cat("Start tuning FDA\n")
    
    if(is.null(ctrl.FDA)){ctrl.FDA <- trControl}
    
    tune.grid <- expand.grid(.degree = 1:2, .nprune = 2:38)
    try(tune.FDA <- caret::train(data@data.env.var, factor(resp), 
                          method = "fda",
                          tuneGrid = tune.grid,                  
                          trControl = ctrl.FDA))
    cat(paste("Finished tuning FDA\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  if('CTA' %in% models){
    cat("Start tuning CTA\n")
    
    if(is.null(ctrl.CTA)){ctrl.CTA <- trControl}    
    
    cat("Tuning Complexity Parameter")    
    try(tune.CTA.rpart <- caret::train(data@data.env.var, resp, 
                                method = "rpart",
                                tuneLength = tuneLength,
                                trControl = ctrl.CTA,
                                metric=metric))
    
    cat("Tuning Max Tree Depth")
    try(tune.CTA.rpart2 <-  caret::train(data@data.env.var, resp,
                                  method = "rpart2",
                                  tuneLength = tuneLength,
                                  trControl = ctrl.CTA,
                                  metric=metric))
    cat(paste("Finished tuning CTA\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  ## Compiling information:
  
  if(!is.null(tune.GLM)){
    tune.GLM
    models.options@GLM$myFormula <- formula(tune.GLM$finalModel)
    models.options@GLM$type <- type.GLM
    models.options@GLM$test <- "no"    
  } else { if('GLM' %in% models){cat("Tuning GLM failed!"); tune.GLM <- "FAILED"}}
  
  if(!is.null(tune.MAXENT.Phillips)){
    if(metric.ME=="ROC"){metric.ME <- "Mean.AUC"}
    if(!metric.ME %in% c("Mean.AUC", "Mean.AUC.DIFF", "delta.AICc")){metric.ME <- "Mean.AUC"; cat("Invalid metric.ME argument! metric.ME was set to Mean.AUC")}
    if(metric.ME == 'Mean.AUC'){
      models.options@MAXENT.Phillips$linear <- grepl("L",tune.MAXENT.Phillips@results[which.max(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
      models.options@MAXENT.Phillips$quadratic <- grepl("Q",tune.MAXENT.Phillips@results[which.max(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
      models.options@MAXENT.Phillips$hinge <- grepl("H",tune.MAXENT.Phillips@results[which.max(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
      models.options@MAXENT.Phillips$product <- grepl("P",tune.MAXENT.Phillips@results[which.max(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
      models.options@MAXENT.Phillips$threshold <- grepl("T",tune.MAXENT.Phillips@results[which.max(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
      models.options@MAXENT.Phillips$beta_threshold <- tune.MAXENT.Phillips@results[which.max(tune.MAXENT.Phillips@results[,metric.ME]),"rm"]  
    }else {       
      models.options@MAXENT.Phillips$linear <- grepl("L",tune.MAXENT.Phillips@results[which.min(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
      models.options@MAXENT.Phillips$quadratic <- grepl("Q",tune.MAXENT.Phillips@results[which.min(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
      models.options@MAXENT.Phillips$hinge <- grepl("H",tune.MAXENT.Phillips@results[which.min(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
      models.options@MAXENT.Phillips$product <- grepl("P",tune.MAXENT.Phillips@results[which.min(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
      models.options@MAXENT.Phillips$threshold <- grepl("T",tune.MAXENT.Phillips@results[which.min(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
      models.options@MAXENT.Phillips$beta_threshold <- tune.MAXENT.Phillips@results[which.min(tune.MAXENT.Phillips@results[,metric.ME]),"rm"]  
    }}else{ if('MAXENT.Phillips' %in% models){cat("Tuning MAXENT.Phillips failed!"); tune.MAXENT.Phillips <- "FAILED"}}
  
  if(!is.null(tune.GAM)){
    if(metric == 'TSS'){
      models.options@GAM$select <- tune.GAM$results[which.max(apply(tune.GAM$results[,c("Sens","Spec")],1,sum)-1),"select"]    
      models.options@GAM$method <- tune.GAM$results[which.max(apply(tune.GAM$results[,c("Sens","Spec")],1,sum)-1),"method"]    
    }else {
      models.options@GAM$select <- tune.GAM$bestTune$select
      models.options@GAM$method <- tune.GAM$bestTune$method
    }}else{ if('GAM' %in% models){cat("Tuning GAM failed!"); tune.GAM <- "FAILED"}}
  
  if(!is.null(tune.GBM)){
    if(metric == 'TSS'){
      models.options@GBM$n.trees <- tune.GBM$results[which.max(apply(tune.GBM$results[,c("Sens","Spec")],1,sum)-1),"n.trees"]
      models.options@GBM$interaction.depth <- tune.GBM$results[which.max(apply(tune.GBM$results[,c("Sens","Spec")],1,sum)-1),"interaction.depth"]    
      models.options@GBM$shrinkage <- tune.GBM$results[which.max(apply(tune.GBM$results[,c("Sens","Spec")],1,sum)-1),"shrinkage"] 
    }else{
      models.options@GBM$n.trees <- tune.GBM$bestTune$n.trees
      models.options@GBM$interaction.depth <- tune.GBM$bestTune$interaction.depth    
      models.options@GBM$shrinkage <- tune.GBM$bestTune$shrinkage 
    }}else{ if('GBM' %in% models){cat("Tuning GBM failed!"); tune.GBM <- "FAILED"}}
  
  if(!is.null(tune.CTA.rpart)){
    if(metric == 'TSS'){
      models.options@CTA$control$cp <- tune.CTA.rpart$results[which.max(apply(tune.CTA.rpart$results[,c("Sens","Spec")],1,sum)-1),"cp"]    
    }else{
      models.options@CTA$control$cp <- tune.CTA.rpart$bestTune      
    }}else{ if('CTA' %in% models){cat("Tuning CTA cp failed!"); tune.CTA.rpart <- "FAILED"}}
  
  if(!is.null(tune.CTA.rpart2)){
    if(metric == 'TSS'){
      models.options@CTA$control$maxdepth <- tune.CTA.rpart2$results[which.max(apply(tune.CTA.rpart2$results[,c("Sens","Spec")],1,sum)-1),"maxdepth"]    
    }else{
      models.options@CTA$control$maxdepth <- tune.CTA.rpart2$bestTune      
    }}else{ if('CTA' %in% models){cat("Tuning CTA maxdepth failed!"); tune.CTA.rpart2 <- "FAILED"}}
  
  if(!is.null(tune.RF)){
    if(metric == 'TSS'){
      models.options@RF$mtry <- tune.RF$results[which.max(apply(tune.RF$results[,c("Sens","Spec")],1,sum)-1),"mtry"]    
    }else{
      models.options@RF$mtry <- tune.RF$bestTune$mtry  
    }}else{ if('RF' %in% models){cat("Tuning RF failed!"); tune.RF <- "FAILED"}}
  
  if(!is.null(tune.ANN)){
    if(metric == 'TSS'){
      models.options@ANN$size <- tune.ANN$results[which.max(apply(tune.ANN$results[,c("Sens","Spec")],1,sum)-1),"size"]    
      models.options@ANN$decay <- tune.ANN$results[which.max(apply(tune.ANN$results[,c("Sens","Spec")],1,sum)-1),"decay"]    
      models.options@ANN$maxit <- maxit.ANN
    }else{
      models.options@ANN$size <- tune.ANN$bestTune$size
      models.options@ANN$decay <- tune.ANN$bestTune$decay 
      models.options@ANN$maxit <- maxit.ANN
    }}else{ if('ANN' %in% models){cat("Tuning ANN failed!"); tune.ANN <- "FAILED"}}
  
  if(!is.null(tune.MARS)){
    if(metric == 'TSS'){
      models.options@ANN$degree <- tune.MARS$results[which.max(apply(tune.MARS$results[,c("Sens","Spec")],1,sum)-1),"degree"]    
      models.options@ANN$nk <- tune.MARS$results[which.max(apply(tune.MARS$results[,c("Sens","Spec")],1,sum)-1),"nk"]    
    }else{
      models.options@MARS$degree <- tune.MARS$bestTune$degree
      models.options@MARS$nk <- tune.MARS$bestTune$nprune
    }}else{ if('MARS' %in% models){cat("Tuning MARS failed!"); tune.MARS <- "FAILED"}}
  
  if(!is.null(tune.FDA)){
    #models.options@FDA$method <- "earth"   
    if(metric == 'TSS'){
      models.options@FDA$add_args <- list(degree=tune.FDA$results[which.max(apply(tune.FDA$results[,c("Sens","Spec")],1,sum)-1),"degree"],   
                                          nk=tune.FDA$results[which.max(apply(tune.FDA$results[,c("Sens","Spec")],1,sum)-1),"nprune"])    
    }else{
      models.options@FDA$add_args <- list(degree=tune.FDA$bestTune$degree,nk=tune.FDA$bestTune$nprune)
    }}else{ if('FDA' %in% models){cat("Tuning FDA failed!"); tune.FDA <- "FAILED"}}
  
  return(list(models.options=models.options,   tune.CTA.rpart = tune.CTA.rpart, tune.CTA.rpart2 = tune.CTA.rpart2,
              tune.RF = tune.RF, tune.ANN = tune.ANN,  tune.MARS = tune.MARS, tune.FDA = tune.FDA, tune.GBM=tune.GBM,
              tune.GAM = tune.GAM, tune.MAXENT.Phillips = tune.MAXENT.Phillips, tune.GLM=tune.GLM))
}


