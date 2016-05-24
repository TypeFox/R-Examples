# Creates a TuneParetoClassifier object encapsulating a classifier.
# <name> is a human-readable name of the classifier.
# <classifier> is the classifier function or its name.
# <classifierParamNames> is a list of possible parameters for <classifier>.
# <predefinedClassifierParams> is a list of default values for parameters of <classifier>.
# <predictor> is a prediction function, or NULL if classifier and predictor are combined.
# <predictorParamNames> is a list of possible parameters for <predictor>.
# <predefinedPredictorParams> is a list of default values for parameters of <predictor>.
# If <useFormula> is true, the classifier is supplied with a formula specifying the relation-
# ship between the data and the classes in a parameter called <formulaName>, and
# the class labels are automatically included in the training data frame.
# <trainDataName>, <trainLabelName>, <testDataName>, and <modelName> specify the names
# of the parameters of <classifier> and <predictor> for the training data, the training labels,
# the test data, and the classification model respectively.
# <requiredPackages> is a character vector of packages to load for the classifier.
tuneParetoClassifier <- function(name, classifier, classifierParamNames=NULL, predefinedClassifierParams=NULL,
                                 predictor=NULL, predictorParamNames=NULL, predefinedPredictorParams=NULL,
                                 useFormula = FALSE, formulaName = "formula",
                                 trainDataName = "x", trainLabelName = "y", testDataName = "newdata",
                                 modelName = "object", requiredPackages=NULL)
{
  # load packages
  if (length(requiredPackages) > 0)
    lapply(requiredPackages,function(package)require(package, character.only=TRUE))
    
  res <- list(name = name,
              classifier = classifier,
              classifierParamNames = classifierParamNames,
              predefinedClassifierParams = predefinedClassifierParams,
              predictor = predictor,
              predictorParamNames = predictorParamNames,
              predefinedPredictorParams = predefinedPredictorParams,
              useFormula = useFormula,
              formulaName = formulaName,
              trainDataName = trainDataName,
              trainLabelName = trainLabelName,
              testDataName = testDataName,
              modelName = modelName,
              requiredPackages = requiredPackages)
  class(res) <- "TuneParetoClassifier"
  return(res)              
}

# TuneParetoClassifier wrapper for CART
tunePareto.tree <- function()
{
  return(tuneParetoClassifier(name = "tree",
                              classifier = "tree", 
                              predictor = "predict", 
                              classifierParamNames = c("weights", "subset",
                                                       "na.action","method","split",
                                                       "mincut", "minsize", "mindev"),
                              predictorParamNames =  c("type"),
                              predefinedPredictorParams = list(type="class"),
                              useFormula=TRUE,
                              formulaName="formula",
                              trainDataName="data",
                              testDataName="newdata",
                              modelName="object",
                              requiredPackages="tree"))
}

# TuneParetoClassifier wrapper for random forests
tunePareto.randomForest <- function()
{
  return(tuneParetoClassifier(name = "randomForest",
                              classifier = "randomForest", 
                              predictor = "predict", 
                              classifierParamNames = c("subset", "na.action",
                                                       "ntree", "mtry",
                                                       "replace", "classwt", 
                                                       "cutoff", "strata",
                                                       "sampsize", "nodesize",
                                                       "maxnodes"),
                              predictorParamNames =  NULL,
                              useFormula=FALSE,
                              trainDataName="x",
                              trainLabelName="y",
                              testDataName="newdata",
                              modelName="object",
                              requiredPackages="randomForest"))
}

# TuneParetoClassifier wrapper for k-Nearest Neighbour
tunePareto.knn <- function()
{
  return(tuneParetoClassifier(name = "knn",
                              classifier = "knn", 
                              predictor = NULL, 
                              classifierParamNames = c("k", "l", "use.all"),
                              predictorParamNames = NULL,
                              useFormula = FALSE,
                              trainDataName = "train",
                              testDataName = "test",
                              trainLabelName = "cl",
                              requiredPackages="class"))
}

# TuneParetoClassifier wrapper for support vector machines
tunePareto.svm <- function()
{
  return(tuneParetoClassifier(name = "svm",
                              classifier = "svm", 
                              predictor = "predict", 
                              classifierParamNames = c("kernel", "degree", "gamma",
                                                       "coef0", "cost", "nu",
                                                       "class.weights", "cachesize", 
                                                       "tolerance", "epsilon",
                                                       "scale", "shrinking", "fitted",
                                                       "subset", "na.action"),
                              useFormula = FALSE,
                              trainDataName = "x",
                              trainLabelName = "y",
                              testDataName = "newdata",
                              modelName = "object",
                              requiredPackages="e1071"))
}

# TuneParetoClassifier wrapper for Naive Bayes classification
tunePareto.NaiveBayes <- function()
{
   return(tuneParetoClassifier(name = "NaiveBayes",
                               classifier = NaiveBayes,
                               predictor  = function(object, newdata, ...)
                                            {
                                              predict(object, newdata, ...)$class
                                            },
                               classifierParamNames =  c("prior","usekernel","fL","subset",
                                                         "na.action","bw","adjust","kernel","weights",
                                                         "window", "width","give.Rkern","n",
                                                         "from","to", "cut","na.rm"),
                               predictorParamNames = "threshold",
                               useFormula = FALSE,
                               trainDataName = "x",
                               testDataName = "newdata",
                               trainLabelName = "grouping",
                               requiredPackages = "klaR"))
} 

# Trains a TuneParetoClassifier object <classifier>
# with the training data supplied in <trainData> and the labels supplied in <trainLabels>.
# Classifier training parameters can be specified in the ... argument.
# Returns a TuneParetoModel object with the training information.
trainTuneParetoClassifier <- function(classifier, trainData, trainLabels, ...)
{
  if(!inherits(classifier, "TuneParetoClassifier"))
    stop("\"classifier\" must be a TuneParetoClassifier object!")
    
  args <- list(...)
  nonmatch <- setdiff(names(args), classifier$classifierParamNames)
  if (length(nonmatch) > 0)
    stop("The following unknown classifier parameters have been specified: ",nonmatch)

  unsetParams <- setdiff(names(classifier$predefinedClassifierParams),names(args))
  args[unsetParams] <- classifier$predefinedClassifierParams[unsetParams]

  trainLabels <- as.factor(trainLabels)
    
  if (is.null(classifier$predictor))
  # empty predictor => combined train/predict method
  {
    res <- list(trainData=trainData, trainLabels=trainLabels, 
                classifierParams = args,
                classifier=classifier)
  }
  else
  # separate training and prediction  
  {
    if (classifier$useFormula)
    {
      trainData <- data.frame(trainLabels,trainData)
      colnames(trainData)[1] <- "Class"
      
      paramList <- list(as.formula("Class ~ ."), trainData)
      names(paramList) <- c(classifier$formulaName, classifier$trainDataName)        
    }
    else
    {  
      paramList <- list(trainData, trainLabels)
      names(paramList) <- c(classifier$trainDataName, classifier$trainLabelName)
    }
    paramList <- c(paramList, args)
  
    # train the classifier
    model <- do.call(classifier$classifier, paramList)
    res <- list(model=model, classifierParams=args, classifier=classifier)
  }
  class(res) <- "TuneParetoModel"
  return(res)
}

# Prediction method for TuneParetoModel objects.
# <object> is a TuneParetoModel training object.
# <newdata> is a matrix or data frame of new samples.
# ... contains additional prediction parameters.
# Returns a vector of predicted class labels.
predict.TuneParetoModel <- function(object, newdata, ...)
{
  if(!inherits(object, "TuneParetoModel"))
    stop("\"object\" must be a TuneParetoModel object!")
    
  args <- list(...)
  nonmatch <- setdiff(names(args), object$classifier$predictorParamNames)
  if (length(nonmatch) > 0)
    stop("The following unknown predictor parameters have been specified: ",nonmatch)
    
  unsetParams <- setdiff(names(object$classifier$predefinedPredictorParams),names(args))
  args[unsetParams] <- object$classifier$predefinedPredictorParams[unsetParams]
  
  if (is.null(object$classifier$predictor))
  # combined classification/prediction method
  {
    
    if (object$classifier$useFormula)
    # build a data frame and a formula
    {
      trainData <- data.frame(object$trainLabels,object$trainData)
      colnames(trainData)[1] <- "Class"
      
      testData <- as.data.frame(newdata)
      colnames(testData) <- colnames(trainData)[-1]
      
      paramList <- list(as.formula("Class ~ ."), trainData, testData)
      names(paramList) <- c(object$classifier$formulaName, 
                            object$classifier$trainDataName, 
                            object$classifier$testDataName)
    }
    else
    # supply the data separately
    {
      paramList <- list(object$trainData, object$trainLabels, newdata)
      names(paramList) <- c(object$classifier$trainDataName, 
                            object$classifier$trainLabelName, 
                            object$classifier$testDataName)
    }
    paramList <- c(paramList, as.list(object$classifierParams))
    
    # call classifier and return predicted labels
    res <- do.call(object$classifier$classifier, paramList)
    if (length(res) != nrow(newdata))
      stop("The length of the prediction vector does not coincide with the number of samples!")
    else
      return(res)  
  }
  else
  # separate training and prediction
  {
    # build parameter list for testing
    paramList <- list(object$model, as.data.frame(newdata))
    names(paramList) <- c(object$classifier$modelName, 
                          object$classifier$testDataName)
    paramList <- c(paramList, args)
    
    # predict the unknown samples
    res <- do.call(object$classifier$predictor, paramList)
    if (length(res) != nrow(newdata))
      stop("The length of the prediction vector does not coincide with the number of samples!")
    else
      return(res)    
  } 
}

# Print method for TuneParetoClassifier objects
print.TuneParetoClassifier <- function(x, ...)
{
  if(!inherits(x, "TuneParetoClassifier"))
    stop("\"x\" must be a TuneParetoClassifier object!")

  cat("TunePareto classifier object\n\n")
  cat("Classifier: ",x$name,"\n",sep="")
  cat("Classifier parameters: ",paste(x$classifierParamNames,collapse=", "),"\n",sep="")
  cat("Predictor parameters: ",paste(x$predictorParamNames,collapse=", "),"\n",sep="")
  cat("Required packages: ",paste(x$requiredPackages,collapse=", "),"\n",sep="")

  return(invisible(x))
}

# Print method for TuneParetoModel objects
print.TuneParetoModel <- function(x, ...)
{
  if(!inherits(x, "TuneParetoModel"))
    stop("\"x\" must be a TuneParetoModel object!")

 cat("TunePareto trained classifier\n\n")
 cat("Classifier: ",x$classifier$name,"\n",sep="") 
 cat("Used parameters: ",
     paste(mapply(function(key,val)
                  {
                    paste(key,val,sep=" = ")
                  },
                  names(x$classifierParams),x$classifierParams),
     collapse=", "),
     "\n",sep="")
 
 return(invisible(x))
}
