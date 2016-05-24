###################################################################################
##                               SelvarLearnLasso.R                              ##
###################################################################################
SelvarLearnLasso <- 
  function(data,  
           knownlabels,
           lambda, 
           rho,
           hybrid.size,  
           models,
           regModel,
           indepModel,
           dataTest,
           labelsTest, 
           nbCores)
  {
    
    # check data parameter
    if(missing(data)){
      stop("data is missing !")
    } 
    if(is.matrix(data) == FALSE & is.data.frame(data) == FALSE) 
      stop(paste(sQuote("data"), "must be a matrix"))
    
    # check lambda parameter
    if(missing(lambda)){
      stop("lambda is missing!")
    }
    if(is.vector(lambda) == FALSE | length(lambda) <= 1){ 
      stop(paste(sQuote("lambda"), "must be a vector with length >= 2"))
    }
    if (sum(lambda<=0)){
      stop("lambda must be greater than 0!")
    }
    
    # check rho parameter
    if(missing(rho)){
      stop("rho is missing!")
    }
    if(is.vector(rho) == FALSE){ 
      stop(paste(sQuote("rho"), "must be a vector"))
    }
    if(sum(rho<=0)){
      stop("rho must be greater than 0!")
    }
    
    
    # check hybrid.size  default value = 3
    if(missing(hybrid.size)){
      hybrid.size <- 3
    }
    if(!is.wholenumber(hybrid.size) | sum(hybrid.size < 1) | hybrid.size > ncol(data)) 
      stop(paste(sQuote("hybrid.size"), "must be a positive integer <= ncol(data)"))
    
    # check models 
    if(missing(models)){
      models <- mixmodGaussianModel(listModels = c("Gaussian_pk_L_C", 
                                                   "Gaussian_pk_Lk_C", 
                                                   "Gaussian_pk_L_Ck", 
                                                   "Gaussian_pk_Lk_Ck"))
    }
    if(!(isS4(models) && is(models, "GaussianModel")))
    {
      stop("models must be an GaussianModel S4 object! (see Rmixmod package)")
    }
    # check criterion parameter
    #     if( sum(criterion %in% c("BIC")) != length(criterion) ){
    #       stop(cat(criterion[which(!(criterion %in% c("ICL")))], "is not a valid criterion name !\n"))
    #     }
    #     
    # check regModel
    if(missing(regModel)){
      regModel <- c("LI", "LB", "LC")
    }
    if( sum(regModel %in% c("LI","LB","LC")) != length(regModel) ){
      stop(cat(regModel[which(!(regModel %in% c("LI","LB","LC")))], "is not a valid regModel name !\n"))
    }
    
    # check indepModel
    if(missing(indepModel)){
      indepModel <- c("LI", "LB")
    }
    if ( sum(indepModel %in% c("LI","LB")) != length(indepModel) ){
      stop(cat(indepModel[which(!(indepModel %in% c("LI","LB")))], "is not a valid indepModel name !\n"))
    }
    
    # check whether the knownLabels is missing
    if( missing(knownlabels)){
      stop("labels are missing!")
    }
    
    if(is.factor(knownlabels))
    {
      knownlabels <- factor(knownlabels, labels = seq(1,length(levels(knownlabels))))
      knownlabels <- as.integer(knownlabels)
    }
    
    # check the number of cluster
    if (min(knownlabels) <= 0 | length(knownlabels) != nrow(data)){
      stop("Each observation in knownLabels must have a valid cluster affectation !")
    }
   
    
    ## check dataTest and labelsTest
    if(missing(dataTest) && !missing(labelsTest))
      stop("dataTest is missing!")
    
    if(!missing(dataTest) && missing(labelsTest))
      stop("labelsTest are missing!")
    
    if(!missing(dataTest))
      if((is.matrix(dataTest) == FALSE & is.data.frame(dataTest) == FALSE)) 
        stop(paste(sQuote("dataTest"), "must be a matrix"))
    
    testing <- TRUE
    if(missing(dataTest) && missing(labelsTest))
      testing <- FALSE
    
    # check nbCores 
    nb.cpus <- detectCores(all.tests = FALSE, logical = FALSE)
    if(missing(nbCores))
    {
      if(nb.cpus > 1)
        nbCores <- 2
      if(nb.cpus == 1)
        nbCores <- 1
    }
    
    
    supervised <- TRUE
    nbCluster <- as.integer(max(knownlabels))
    criterion <- "BIC" 
    data <- as.matrix(data)
    n <- as.integer(nrow(data))
    p <- as.integer(ncol(data))
    OrderVariable <- rep(NA, p) 
    dataStand <- scale(data, TRUE, TRUE)
    print("............... start  variable  ranking .................................... ")
    OrderVariable <- SortvarLearn(dataStand, knownlabels, lambda, rho, nbCores)
    print("................. variable ranking .... done ................................ ") 
    bestModel <- list()
    print(" ...... SRUW  selection with BIC criterion ...... ")
    VariableSelectRes <- VariableSelection(data,
                                           nbCluster,
                                           models,
                                           criterion,
                                           OrderVariable,
                                           hybrid.size,
                                           supervised,
                                           knownlabels,
                                           nbCores)
    print(" ...... model selection with BIC criterion......")
    bestModel <- ModelSelectionClust(VariableSelectRes,
                                     data,
                                     regModel,
                                     indepModel,
                                     nbCores)
    
    if(testing)
    {
      dataAux <- as.data.frame(data[,bestModel$S])
      dataTestAux <- as.data.frame(dataTest[,bestModel$S])
      model <- mixmodGaussianModel(listModels = bestModel$model)
      learn <- mixmodLearn(dataAux, knownLabels = knownlabels, models = model)
      predict <- mixmodPredict(dataTestAux, classificationRule = learn["bestResult"])
      bestModel$proba <- predict["proba"]
      bestModel$partition <- predict["partition"]
      bestModel$error <- 1 - mean(predict["partition"] == labelsTest)
    }
    else
    {
      bestModel$proba <- NULL
      bestModel$partition <- NULL
      bestModel$error <- NULL
    }
    
    
    return(bestModel)  
    
  }