#################################################################
# Objective functions and utilities for objective functions
#################################################################

# Creates a new objective from the specified objective
# function <objectiveFunction> and its precalculation
# method <precalculationFunction> with additional parameters
# <precalculationParams>. 
# <direction> specifies whether the objective is 
# minimized or maximized. 
# <name> is a short description of the objective.
#
# Returns an object of class TuneParetoObjective
createObjective <- function(precalculationFunction,
                            precalculationParams = NULL,
                            objectiveFunction,
                            objectiveFunctionParams = NULL,
                            direction=c("minimize","maximize"),
                            name)
{
  res <- list(precalculationFunction = match.fun(precalculationFunction),
              precalculationParams = precalculationParams,
              objectiveFunction = match.fun(objectiveFunction),
              objectiveFunctionParams = objectiveFunctionParams,
              minimize = (match.arg(direction, c("minimize","maximize")) == "minimize"),
              name = name)
  class(res) <- "TuneParetoObjective"
  return(res)
}

# Groups all objectives in <objectiveFunctionList> that
# use the same preprocessing with the same parameters.
# This speeds up calculations.
groupByPrecalculation <- function(objectiveFunctionList)
{
  # determine unique set of precalculation functions
  funcs  <- sapply(objectiveFunctionList,function(obj)obj$precalculationFunction)
  uniqueFuncs <- unique(funcs)
  resultList <- list()
  indices <- c()
  
  for (func in uniqueFuncs)
  {
    # find objectives with the same preprocessing as <func>
    newIndices <- which(sapply(funcs, function(f)identical(f,func)))
    indices <- c(indices,newIndices)
    
    # find unique set of parameters for the preprocessing
    subset <- objectiveFunctionList[newIndices]
    params <- lapply(subset,function(obj)obj$precalculationParams)
    uniqueParams <- unique(params)
    newObjectives <- lapply(uniqueParams,function(param)
                            {
                              # join objectives with the same parameters
                              toJoin <- subset[sapply(params, function(p)identical(p,param))]
                              
                              res <- toJoin[[1]]

                              if (length(toJoin) > 1)
                              {
                                res$objectiveFunctionParams[[1]] <- list(res$objectiveFunctionParams)
                                for (i in 2:length(toJoin))
                                {
                                  res$objectiveFunction <- c(res$objectiveFunction, toJoin[[i]]$objectiveFunction)
                                  res$objectiveFunctionParams[[i]] <- 
                                        list(toJoin[[i]]$objectiveFunctionParams)
                                  res$minimize <- c(res$minimize, toJoin[[i]]$minimize)
                                  res$name <- c(res$name, toJoin[[i]]$name)
                                }
                              }
                              res
                            })
                            
    # add new joint objective to result list
    resultList[(length(resultList) + 1):(length(resultList) + length(newObjectives))] <- newObjectives
  }
  
  # calculate permutation of original indices for subsequent reordering
  permutation <- indices
  for (i in seq_along(permutation))
    permutation[indices[i]] <- i
    
  return(list(grouping=resultList, permutation=permutation))
  
}

# Predefined precalculation function that trains the classifier on the whole
# data set and predicts the same samples with the trained classifier.
# For parameters, see callClassifier.
#
# Returns a list containing a vector of predicted labels and a vector of true labels
reclassification <- function(data, labels, classifier, classifierParams, predictorParams)
{
  arglist <- list(classifier=classifier,trainData=data,trainLabels=labels)
  arglist <- c(arglist, classifierParams)
  train <- do.call(trainTuneParetoClassifier,arglist)
  
  arglist <- list(object=train, newdata=data)
  arglist <- c(arglist, predictorParams)
  
  predictedLabs <- do.call(predict,arglist)
  
  #labels <- as.integer(as.character(labels))

  res <- list(predictedLabels=predictedLabs, trueLabels=labels, model=train)
  class(res) <- "ClassificationResult"
                             
  return(res)
}

# Generate a list of partitions for a cross-validation that can
# be used for all cross-validation methods. 
# For parameters, see crossValidation
generateCVRuns <- function(labels, ntimes = 10, nfold = 10, 
                           leaveOneOut=FALSE, stratified = FALSE)
{
  if(leaveOneOut)
	ntimes <- 1

  numSamples <- length(labels)
  res <- lapply(1:ntimes,function(run)
  # for each run
  {
    # calculate folds
	  if (leaveOneOut)
	  {
		  indices <- as.list(1:numSamples)
	  }
	  else
	  {
		  if(stratified)
		  {
			  classes <- unique(labels)
			  sing.perm <- lapply(classes, function(cl){
				  index <- which(labels == cl)
				  sample(index, length(index))
			  })
			  permut <- unlist(sing.perm)
			  indices <- lapply(1:nfold,function(i){c()})
			  for(i in 1:numSamples)
			  {
				  k = i%%nfold
				  if(k==0)
				   k = nfold
				
				  indices[[k]] <- c(indices[[k]], permut[i])
			  }
		  }
		  else
		  {
		    # permute the indices of the samples
		    permut <- sample(1:numSamples, numSamples,replace=FALSE)
		    indices <- lapply(1:nfold, function(i)
		    {
			    # split the samples in nfold groups
			    permut[seq(i, numSamples, nfold)]
		    })
		  }
	  }
	  names(indices) <- paste("Fold ",1:nfold)
	  return(indices)
	})
	names(res) <- paste("Run ", 1:ntimes)
	return(res)
}

# Predefined precalculation function that performs a cross-validation on <data>.
# <ntimes> is the number of repetitions of the cross-validation.
# <nfold> is the number of groups in each cross-validation run.
# If <leaveOneOut> is true, a leave-one-out cross-validation is performed.
# If <stratified> is true, a stratified cross-validation is carried out.
# For further parameters, see tuneParetoClassifier.
#
# Returns a list containing a sub-list for each run. Each of these sub-lists contains
# a vector of true labels and predicted labels for each fold.
crossValidation <- function(data, labels, classifier, classifierParams, predictorParams,
                            ntimes = 10, nfold = 10, leaveOneOut=FALSE, stratified = FALSE,
                            foldList = NULL)
{
    
  if (is.null(foldList))
    foldList <- generateCVRuns(labels=labels, 
                            ntimes=ntimes, nfold=nfold,
                            leaveOneOut=leaveOneOut,
                            stratified=stratified)
  
  res <- lapply(foldList,function(indices)
  # for each run
  {    
	  
	  return(lapply(indices, function(fold)
	  # for each fold
	  {
	    # split up data
	    trainData <- data[-fold,,drop=FALSE]
	    trainLabels <- labels[-fold]
	    testData <- data[fold,,drop=FALSE]
	    testLabels <- labels[fold]
	    
	    arglist <- list(classifier=classifier,trainData=trainData,trainLabels=trainLabels)
      arglist <- c(arglist, classifierParams)
      train <- do.call(trainTuneParetoClassifier,arglist)
  
      arglist <- list(object=train, newdata=testData)
      arglist <- c(arglist, predictorParams)
  
	    # predict test data
      res1 <- do.call(predict,arglist)
	        
      # return the true labels and the predicted labels
      res1 <- list(predictedLabels=res1, trueLabels=testLabels, model=train)
      class(res1) <- "ClassificationResult"
	    return(res1)
	  }))
	})
	return(res)
}

# Predefined objective calculating the accuracy
# of a reclassification experiment
reclassAccuracy <- function(saveModel = FALSE)
{
  createObjective(precalculationFunction = reclassification,
                  precalculationParams = NULL,
                  objectiveFunction = function(result, saveModel)
                  {
                    res <- sum(result$predictedLabels == result$trueLabels, na.rm = TRUE)/length(result$trueLabels)
                    if (saveModel)
                    {
                      return(list(additionalData = result$model, fitness = res))
                    }
                    else
                      return(res)                                            
                  },
                  objectiveFunctionParams = list(saveModel = saveModel),
                  direction = "maximize",
                  name = "Reclass.Accuracy")
}

# Predefined objective calculating the error percentage
# of a reclassification experiment
reclassError <- function(saveModel=FALSE)
{
  createObjective(precalculationFunction = reclassification,
                  precalculationParams = NULL,
                  objectiveFunction = function(result, saveModel)
                  {
                    res <- sum(is.na(result$predictedLabels) | 
                                 result$predictedLabels != result$trueLabels) / 
                      length(result$trueLabels)
                    if (saveModel)
                    {
                      return(list(additionalData=result$model, fitness=res))
                    }
                    else
                      return(res)                                             
                  },
                  objectiveFunctionParams = list(saveModel=saveModel),                                      
                  direction="minimize",
                  name="Reclass.Error")
}

# Predefined objective calculating the weighted error percentage
# of a reclassification experiment 
reclassWeightedError <- function(saveModel=FALSE)
{
  createObjective(precalculationFunction = reclassification,
                  precalculationParams = NULL,
                  objectiveFunction = function(result, saveModel)
                  {
                    classes <- sort(unique(result$trueLabels))
                    res <- sum(sapply(classes,function(cl)
                    {
                      sum(result$trueLabels == cl &
                            result$predictedLabels != cl)/
                        sum(result$trueLabels == cl)
                    })) / length(classes)
                    if (saveModel)
                    {
                      return(list(additionalData=result$model, fitness=res))
                    }
                    else
                      return(res)       
                  },
                  objectiveFunctionParams = list(saveModel=saveModel),                                      
                  direction="minimize",
                  name="Reclass.WeightedError")
}

# Predefined objective calculating the sensitivity
# of a reclassification experiment
reclassSensitivity <- function(caseClass, saveModel=FALSE)
{
  createObjective(precalculationFunction = reclassification,
                  precalculationParams = NULL,
                  objectiveFunction = function(result, caseClass, saveModel)
                                      {
                                        res <- sum(result$predictedLabels[result$trueLabels == caseClass] 
                                                    == caseClass,na.rm=TRUE) / 
                                               sum(result$trueLabels == caseClass)
                                        if (saveModel)
                                        {
                                          return(list(additionalData=result$model, fitness=res))
                                        }
                                        else
                                          return(res)                
                                      },
                  objectiveFunctionParams = list(caseClass=caseClass, saveModel=saveModel),
                  direction="maximize",
                  name="Reclass.Sensitivity")
}

# Predefined objective calculating the recall
# of a reclassification experiment, which is the same as
# the sensitivity of the experiment
reclassRecall <- function(caseClass, saveModel = FALSE)
{
  objective <- reclassSensitivity(caseClass = caseClass, saveModel = saveModel)
  objective$name <- "Reclass.Recall"
  return(objective)
}

# Predefined objective calculating the true positive rate
# of a reclassification experiment, which is the same as
# the sensitivity of the experiment
reclassTruePositive <- function(caseClass, saveModel = FALSE)
{
  objective <- reclassSensitivity(caseClass = caseClass, saveModel = saveModel)
  objective$name <- "Reclass.TruePositive"
  return(objective)
}

# Predefined objective calculating the specificity
# of a reclassification experiment
reclassSpecificity <- function(caseClass, saveModel=FALSE)
{
  createObjective(precalculationFunction = reclassification,
                  precalculationParams = NULL,
                  objectiveFunction = function(result, caseClass, saveModel)
                                      {
                                        res <- sum(result$predictedLabels[result$trueLabels != caseClass] 
                                                    != caseClass,na.rm=TRUE) / 
                                               sum(result$trueLabels != caseClass)
                                        if (saveModel)
                                        {
                                          return(list(additionalData=result$model, fitness=res))
                                        }
                                        else
                                          return(res)                                               
                                      },
                  objectiveFunctionParams = list(caseClass=caseClass, saveModel=saveModel),                                      
                  direction="maximize",
                  name="Reclass.Specificity")
}

# Predefined objective calculating the true negative rate
# of a reclassification experiment, which is the same as
# the specificity of the experiment
reclassTrueNegative <- function(caseClass, saveModel = FALSE)
{
  objective <- reclassSpecificity(caseClass = caseClass, saveModel = saveModel)
  objective$name <- "Reclass.TrueNegative"
  return(objective)
}

# Predefined objective calculating the fallout
# of a reclassification experiment
reclassFallout <- function(caseClass, saveModel = FALSE)
{
  createObjective(precalculationFunction = reclassification,
                  precalculationParams = NULL,
                  objectiveFunction = function(result, caseClass, saveModel)
                  {
                    case <- result$predictedLabels == caseClass
                    case[is.na(case)] <- TRUE
                    res <- sum(case[result$trueLabels != caseClass])/sum(result$trueLabels != caseClass)
                    if (saveModel)
                    {
                      return(list(additionalData = result$model, fitness = res))
                    }
                    else
                      return(res)                                           
                  },
                  objectiveFunctionParams = list(caseClass = caseClass, saveModel = saveModel),
                  direction = "minimize",
                  name = "Reclass.Fallout")
}

# Predefined objective calculating the false positive rate
# of a reclassification experiment, which is the same as
# the fallout of the experiment
reclassFalsePositive <- function(caseClass, saveModel = FALSE) 
{
  objective <- reclassFallout(caseClass = caseClass, saveModel = saveModel)
  objective$name <- "Reclass.FalsePositive"
  return(objective)
}

# Predefined objective calculating the miss rate
# of a reclassification experiment
reclassMiss <- function(caseClass, saveModel = FALSE)
{
  createObjective(precalculationFunction = reclassification,
                  precalculationParams = NULL,
                  objectiveFunction = function(result, caseClass, saveModel)
                  {
                    nCase <- result$predictedLabels != caseClass
                    nCase[is.na(nCase)] <- TRUE
                    res <- sum(nCase & result$trueLabels == caseClass)/sum(result$trueLabels == caseClass)
                    if (saveModel)
                    {
                      return(list(additionalData = result$model, fitness = res))
                    }
                    else
                      return(res)                                            
                  },
                  objectiveFunctionParams = list(caseClass = caseClass, saveModel = saveModel),
                  direction = "minimize",
                  name = "Reclass.Miss")
}

# Predefined objective calculating the false negative rate
# of a reclassification experiment, which is the same as
# the miss rate of the experiment
reclassFalseNegative <- function(caseClass, saveModel = FALSE) 
{
  objective <- reclassMiss(caseClass = caseClass, saveModel = saveModel)
  objective$name <- "Reclass.FalseNegative"
  return(objective)
}

# Predefined objective calculating the precision
# of a reclassification experiment
reclassPrecision <- function(caseClass, saveModel = FALSE)
{
  createObjective(precalculationFunction = reclassification,
                  precalculationParams = NULL,
                  objectiveFunction = function(result, caseClass, saveModel)
                  {
                    case <- result$predictedLabels == caseClass
                    case[is.na(case) & !case] <- TRUE
                    res <- sum(result$predictedLabels == result$trueLabels & result$trueLabels == caseClass, na.rm = TRUE)/sum(case, na.rm = TRUE)
                    if (saveModel)
                    {
                      return(list(additionalData = result$model, fitness = res))
                    }
                    else
                      return(res)
                  },
                  objectiveFunctionParams = list(caseClass = caseClass, saveModel = saveModel),
                  direction = "maximize",
                  name = "Reclass.Precision")
}

# Predefined objective calculating the positive predictive value
# of a reclassification experiment, which is the same as
# the precision of the experiment
reclassPPV <- function(caseClass, saveModel = FALSE)
{
  objective <- reclassPrecision(caseClass = caseClass, saveModel = saveModel)
  objective$name <- "Reclass.PPV"
  return(objective)
}

# Predefined objective calculating the negative predictive value
# of a reclassification experiment
reclassNPV <- function(caseClass, saveModel = FALSE)
{
  createObjective(precalculationFunction = reclassification,
                  precalculationParams = NULL,
                  objectiveFunction = function(result, caseClass, saveModel)
                  {
                    nCase <- result$predictedLabels != caseClass
                    nCase[is.na(nCase)] <- TRUE
                    res <- sum(result$predictedLabels != caseClass & result$trueLabels != caseClass, na.rm = TRUE)/ sum(nCase)
                    if (saveModel)
                    {
                      return(list(additionalData = result$model, fitness = res))
                    }
                    else
                      return(res)                                            
                  },
                  objectiveFunctionParams = list(caseClass = caseClass, saveModel = saveModel),
                  direction = "maximize",
                  name = "Reclass.NPV")
}

# Predefined objective calculating the confusion of two classes
# in a reclassification experiment
reclassConfusion <- function(trueClass, predictedClass, saveModel=FALSE)
{
  createObjective(precalculationFunction = reclassification,
                  precalculationParams = NULL,
                  objectiveFunction = function(result, trueClass, predictedClass, saveModel)
                                      {
                                        res <- sum(result$predictedLabels[result$trueLabels == trueClass]
                                                    == predictedClass,na.rm=TRUE) / 
                                               sum(result$trueLabels == trueClass)
                                        if (saveModel)
                                        {
                                          return(list(additionalData=result$model, fitness=res))
                                        }
                                        else
                                          return(res)                                               
                                      },
                  objectiveFunctionParams = list(trueClass=trueClass, predictedClass=predictedClass, saveModel=saveModel),
                  direction="minimize",
                  name="Reclass.Confusion")
}

# Predefined objective calculating the accuracy
# of a cross-validation experiment
cvAccuracy <- function (nfold = 10, ntimes = 10, leaveOneOut = FALSE, stratified = FALSE, foldList = NULL, saveModel = FALSE) 
  
{
  createObjective(precalculationFunction = crossValidation,
                  precalculationParams = list(nfold = nfold, ntimes = ntimes, leaveOneOut = leaveOneOut, stratified = stratified, foldList = foldList),
                  objectiveFunction = function(result, saveModel)
                  {
                    numSamples <- sum(sapply(result[[1]], function(fold) length(fold$trueLabels)))
                    res <- mean(sapply(result, function(run)
                    {
                      predictedLabels <- unlist(lapply(run, function(fold) fold$predictedLabels))
                      trueLabels <- unlist(lapply(run, function(fold) fold$trueLabels))
                      return(sum(predictedLabels == trueLabels, na.rm = TRUE) / numSamples)
                    }))
                    if (saveModel)
                    {
                      model <- lapply(result, function(run)
                      {
                        lapply(run, function(fold)
                        {
                          fold$model
                        })
                      })
                      return(list(additionalData = model, fitness = res))
                    }
                    else
                      return(res)                                            
                  },
                  objectiveFunctionParams = list(saveModel = saveModel),
                  direction = "maximize",
                  name = "CV.Accuracy")
}

# Predefined objective calculating the error percentage
# of a cross-validation experiment
cvError <- function(nfold=10, ntimes=10, leaveOneOut=FALSE, stratified = FALSE, foldList=NULL, saveModel=FALSE)
{
    createObjective(precalculationFunction = crossValidation,
                    precalculationParams = list(nfold = nfold, 
                                                ntimes = ntimes, 
                                                leaveOneOut = leaveOneOut,
                                                stratified = stratified,
                                                foldList = foldList),
                    objectiveFunction = function(result, saveModel)
                                        {
                                          numSamples <- sum(sapply(result[[1]], function(fold)length(fold$trueLabels)))
                                          res <- mean(sapply(result,function(run)
                                                      {
                                                        sum(unlist(lapply(run,function(fold)
                                                        {
                                                          is.na(fold$predictedLabels) | 
                                                          fold$predictedLabels != fold$trueLabels
                                                        })))
                                                      }))/numSamples
                                          if (saveModel)
                                          {
                                            model <- lapply(result,function(run)
                                                      {
                                                        lapply(run,function(fold)
                                                        {
                                                          fold$model
                                                        })
                                                      })
                                            return(list(additionalData=model, fitness=res))
                                          }
                                          else
                                            return(res)                                                   
                                          },
                    objectiveFunctionParams = list(saveModel=saveModel),                                        
                    direction="minimize",
                    name="CV.Error")
}

# Predefined objective calculating the variance of the error percentage
# in a cross-validation experiment
cvErrorVariance <- function(nfold=10, ntimes=10, leaveOneOut=FALSE, stratified = FALSE, foldList=NULL, saveModel=FALSE)
{
    createObjective(precalculationFunction = crossValidation,
                    precalculationParams = list(nfold = nfold, 
                                                ntimes = ntimes, 
                                                leaveOneOut = leaveOneOut,
                                                stratified = stratified,
                                                foldList = foldList),
                    objectiveFunction = function(result, saveModel)
                                        {
                                          numSamples <- sum(sapply(result[[1]], function(fold)length(fold$trueLabels)))
                                          res <- var(sapply(result,function(run)
                                                      {
                                                        sum(unlist(lapply(run,function(fold)
                                                        {
                                                          is.na(fold$predictedLabels) | 
                                                          fold$predictedLabels != fold$trueLabels
                                                        })))
                                                      }))/(numSamples^2)
                                          if (saveModel)
                                          {
                                            model <- lapply(result,function(run)
                                                      {
                                                        lapply(run,function(fold)
                                                        {
                                                          fold$model
                                                        })
                                                      })
                                            return(list(additionalData=model, fitness=res))
                                          }
                                          else
                                            return(res)                                                      
                                        },
                    objectiveFunctionParams = list(saveModel=saveModel),                                        
                    direction="minimize",
                    name="CV.Error.Var")
}

# Predefined objective calculating the error percentage
# of a cross-validation experiment
cvWeightedError <- function(nfold=10, ntimes=10, leaveOneOut=FALSE, stratified = FALSE, foldList=NULL, saveModel=FALSE)
{
    createObjective(precalculationFunction = crossValidation,
                    precalculationParams = list(nfold = nfold, 
                                                ntimes = ntimes, 
                                                leaveOneOut = leaveOneOut,
                                                stratified = stratified,
                                                foldList = foldList),
                    objectiveFunction = function(result, saveModel)
                                        {
                                          res <- mean(sapply(result,
                                                 function(run)
                                          {
                                            predictedLabels <- unlist(lapply(run,
                                                                      function(fold)fold$predictedLabels))
                                            trueLabels <- unlist(lapply(run,
                                                                 function(fold)fold$trueLabels))
                                            classes <- sort(unique(trueLabels))
                                            
                                            sum(sapply(classes,function(cl)
                                              {
                                                sum(trueLabels == cl &
                                                    predictedLabels != cl) /
                                                sum(trueLabels == cl)                                                    
                                                
                                              })) / length(classes)
                                          }))                                          
                                          if (saveModel)
                                          {
                                            model <- lapply(result,function(run)
                                                      {
                                                        lapply(run,function(fold)
                                                        {
                                                          fold$model
                                                        })
                                                      })
                                            return(list(additionalData=model, fitness=res))
                                          }
                                          else
                                            return(res)
                                        },
                    objectiveFunctionParams = list(saveModel=saveModel),                                        
                    direction="minimize",
                    name="CV.WeightedError")
}

# Predefined objective calculating the sensitivity
# of a cross-validation experiment
cvSensitivity <- function(nfold=10, ntimes=10, leaveOneOut=FALSE, stratified = FALSE, foldList=NULL, caseClass, saveModel=FALSE)
{
  createObjective(precalculationFunction = crossValidation,
                  precalculationParams = list(nfold = nfold, 
                                              ntimes = ntimes, 
                                              leaveOneOut = leaveOneOut,
                                              stratified = stratified,
                                              foldList = foldList),
                  objectiveFunction = function(result, caseClass, saveModel)
                  {
                    res <- mean(sapply(result,
                                       function(run)
                                       {
                                         predictedLabels <- unlist(lapply(run,
                                                                          function(fold)fold$predictedLabels))
                                         trueLabels <- unlist(lapply(run,
                                                                     function(fold)fold$trueLabels))
                                         
                                         return(sum(predictedLabels[trueLabels == caseClass] 
                                                    == caseClass,na.rm=TRUE) / 
                                                  sum(trueLabels == caseClass))
                                       }))
                    if (saveModel)
                    {
                      model <- lapply(result,function(run)
                      {
                        lapply(run,function(fold)
                        {
                          fold$model
                        })
                      })
                      return(list(additionalData=model, fitness=res))
                    }
                    else
                      return(res)
                  },
                  objectiveFunctionParams = list(caseClass=caseClass, saveModel=saveModel),                                        
                  direction="maximize",
                  name="CV.Sensitivity")
}

# Predefined objective calculating the recall
# of a cross-validation experiment, which is the same as
# the sensitivity of the experiment
cvRecall <- function(nfold = 10, ntimes = 10, leaveOneOut = FALSE, stratified = FALSE, foldList = NULL, caseClass, saveModel = FALSE)
{
  objective <- cvSensitivity(nfold = nfold, ntimes = ntimes, leaveOneOut = leaveOneOut, stratified = stratified, foldList = foldList, caseClass = caseClass, saveModel = saveModel)
  objective$name <- "CV.Recall"
  return(objective)
}

# Predefined objective calculating the true positive rate
# of a cross-validation experiment, which is the same as
# the sensitivity of the experiment
cvTruePositive <- function (nfold = 10, ntimes = 10, leaveOneOut = FALSE, stratified = FALSE, foldList = NULL, caseClass, saveModel = FALSE)
{
  objective <- cvSensitivity(nfold = nfold, ntimes = ntimes, leaveOneOut = leaveOneOut, stratified = stratified, foldList = foldList, caseClass = caseClass, saveModel = saveModel)
  objective$name <- "CV.TruePositive"
  return(objective)
}

# Predefined objective calculating the specificity
# of a cross-validation experiment
cvSpecificity <- function(nfold=10, ntimes=10, leaveOneOut=FALSE, stratified = FALSE, foldList=NULL, caseClass, saveModel=FALSE)
{
  createObjective(precalculationFunction = crossValidation,
                  precalculationParams = list(nfold = nfold, 
                                              ntimes = ntimes, 
                                              leaveOneOut = leaveOneOut,
                                              stratified = stratified,
                                              foldList = foldList),
                  objectiveFunction = function(result, caseClass, saveModel)
                  {
                    res <- mean(sapply(result,
                                       function(run)
                                       {
                                         predictedLabels <- unlist(lapply(run,
                                                                          function(fold)fold$predictedLabels))
                                         trueLabels <- unlist(lapply(run,
                                                                     function(fold)fold$trueLabels))
                                         return(sum(predictedLabels[trueLabels != caseClass] 
                                                    != caseClass,na.rm=TRUE) / 
                                                  sum(trueLabels != caseClass))
                                       }))
                    if (saveModel)
                    {
                      model <- lapply(result,function(run)
                      {
                        lapply(run,function(fold)
                        {
                          fold$model
                        })
                      })
                      return(list(additionalData=model, fitness=res))
                    }
                    else
                      return(res)
                  },
                  objectiveFunctionParams = list(caseClass=caseClass, saveModel=saveModel),                                        
                  direction="maximize",
                  name="CV.Specificity")
}

# Predefined objective calculating the true negative rate
# of a cross-validation experiment, which is the same as
# the specificity of the experiment
cvTrueNegative <- function (nfold = 10, ntimes = 10, leaveOneOut = FALSE, stratified = FALSE, foldList = NULL, caseClass, saveModel = FALSE)
{
  objective <- cvSpecificity(nfold = nfold, ntimes = ntimes, leaveOneOut = leaveOneOut, stratified = stratified, foldList = foldList, caseClass = caseClass, saveModel = saveModel)
  objective$name <- "CV.TrueNegative"
  return(objective)
}

# Predefined objective calculating the fallout
# of a cross-validation experiment
cvFallout <- function(nfold = 10, ntimes = 10, leaveOneOut = FALSE, stratified = FALSE, foldList = NULL, caseClass, saveModel = FALSE)
{
  createObjective(precalculationFunction = crossValidation,
                  precalculationParams = list(nfold = nfold, ntimes = ntimes, leaveOneOut = leaveOneOut, stratified = stratified, foldList = foldList),
                  objectiveFunction = function(result, caseClass, saveModel)
                  {
                    res <-mean(sapply(result, function(run)
                    {
                      predictedLabels <- unlist(lapply(run, function(fold) fold$predictedLabels))
                      trueLabels <- unlist(lapply(run, function(fold) fold$trueLabels))                                               
                      case <- predictedLabels == caseClass
                      case[is.na(case)] <- TRUE
                      return(sum(case[trueLabels != caseClass])/sum(trueLabels != caseClass))
                    }))
                    if (saveModel)
                    {
                      model <- lapply(result, function(run)
                      {
                        lapply(run, function(fold)
                        {
                          fold$model
                        })
                      })
                      return(list(additionalData = model, fitness = res))
                    }
                    else
                      return(res)                                            
                  },
                  objectiveFunctionParams = list(caseClass = caseClass, saveModel = saveModel),
                  direction = "minimize",
                  name = "CV.Fallout")
}

# Predefined objective calculating the false positive rate
# of a cross-validation experiment, which is the same as
# the fallout of the experiment
cvFalsePositive <- function(nfold = 10, ntimes = 10, leaveOneOut = FALSE, stratified = FALSE, foldList = NULL, caseClass, saveModel = FALSE)
{
  objective <- cvFallout(nfold = 10, ntimes = 10, leaveOneOut = FALSE, stratified = FALSE, foldList = NULL, caseClass, saveModel = saveModel)
  objective$name <- "CV.FalsePositive"
  return(objective)
}

# Predefined objective calculating the miss rate
# of a cross-validation experiment
cvMiss <- function(nfold = 10, ntimes = 10, leaveOneOut = FALSE, stratified = FALSE, foldList = NULL, caseClass, saveModel = FALSE)
{
  createObjective(precalculationFunction = crossValidation,
                  precalculationParams = list(nfold = nfold, ntimes = ntimes, leaveOneOut = leaveOneOut, stratified = stratified, foldList = foldList),
                  objectiveFunction = function(result, caseClass, saveModel)
                  {
                    res <- mean(sapply(result, function(run)
                    {
                      predictedLabels <- unlist(lapply(run, function(fold) fold$predictedLabels))
                      trueLabels <- unlist(lapply(run, function(fold) fold$trueLabels))
                      nCase <- predictedLabels != caseClass
                      nCase[is.na(nCase)] <- TRUE
                      return(sum(nCase & trueLabels == caseClass)/sum(trueLabels == caseClass))
                    }))
                    if (saveModel)
                    {
                      model <- lapply(result, function(run)
                      {
                        lapply(run, function(fold)
                        {
                          fold$model
                        })
                      })
                      return(list(additionalData = model, fitness = res))
                    }
                    else
                      return(res)                                            
                  },
                  objectiveFunctionParams = list(caseClass = caseClass, saveModel = saveModel),
                  direction = "minimize",
                  name = "CV.Miss")
}

# Predefined objective calculating the false negative rate
# of a cross-validation experiment, which is the same as
# the miss rate of the experiment
cvFalseNegative <- function(nfold = 10, ntimes = 10, leaveOneOut = FALSE, stratified = FALSE, foldList = NULL, caseClass, saveModel = FALSE)
{
  objective <- cvMiss(nfold = 10, ntimes = 10, leaveOneOut = FALSE, stratified = FALSE, foldList = NULL, caseClass = caseClass, saveModel = saveModel)
  objective$name <- "CV.FalseNegative"
  return(objective)
}

# Predefined objective calculating the precision
# of a cross-validation experiment
cvPrecision <- function(nfold = 10, ntimes = 10, leaveOneOut = FALSE, stratified = FALSE, foldList = NULL, caseClass, saveModel = FALSE)
{
  createObjective(precalculationFunction = crossValidation,
                  precalculationParams = list(nfold = nfold, ntimes = ntimes, leaveOneOut = leaveOneOut, stratified = stratified, foldList = foldList),
                  objectiveFunction = function(result, caseClass, saveModel)
                  {
                    res <- mean(sapply(result, function(run)
                    {
                      predictedLabels <- unlist(lapply(run, function(fold) fold$predictedLabels))
                      trueLabels <- unlist(lapply(run, function(fold) fold$trueLabels))
                      case <- predictedLabels == caseClass
                      case[is.na(case) & !case] <- TRUE
                      return(sum(predictedLabels == trueLabels & trueLabels == caseClass, na.rm = TRUE)/sum(case, na.rm = TRUE))
                    }))
                    if (saveModel)
                    {
                      model <- lapply(result, function(run)
                      {
                        lapply(run, function(fold)
                        {
                          fold$model
                        })
                      })
                      return(list(additionalData = model, fitness = res))
                    }
                    else
                      return(res)                                            
                  },
                  objectiveFunctionParams = list(caseClass = caseClass, saveModel = saveModel),
                  direction = "maximize",
                  name = "CV.Precision")
}

# Predefined objective calculating the positive predictive value
# of a cross-validation experiment, which is the same as
# the precision of the experiment
cvPPV <- function(nfold = 10, ntimes = 10, leaveOneOut = FALSE, stratified = FALSE, foldList = NULL, caseClass, saveModel = FALSE)
{
  objective <- cvPrecision(nfold = nfold, ntimes = ntimes, leaveOneOut = leaveOneOut, stratified = stratified, foldList = foldList, caseClass, saveModel = saveModel)
  objective$name <- "CV.PPV"
  return(objective)
}

# Predefined objective calculating the negative predictive value
# of a cross-validation experiment
cvNPV <- function(nfold = 10, ntimes = 10, leaveOneOut = FALSE, stratified = FALSE, foldList = NULL, caseClass, saveModel = FALSE)
{
  createObjective(precalculationFunction = crossValidation,
                  precalculationParams = list(nfold = nfold, ntimes = ntimes, leaveOneOut = leaveOneOut, stratified = stratified, foldList = foldList),
                  objectiveFunction = function(result, caseClass, saveModel)
                  {
                    res <- mean(sapply(result, function(run)
                    {
                      predictedLabels <- unlist(lapply(run, function(fold) fold$predictedLabels))
                      trueLabels <- unlist(lapply(run, function(fold) fold$trueLabels))
                      nCase <- predictedLabels != caseClass
                      nCase[is.na(nCase)] <- TRUE            
                      return(sum(predictedLabels != caseClass & trueLabels != caseClass, na.rm = TRUE)/ sum(nCase))
                    }))
                    if (saveModel)
                    {
                      model <- lapply(result, function(run)
                      {
                        lapply(run, function(fold)
                        {
                          fold$model
                        })
                      })
                      return(list(additionalData = model, fitness = res))
                    }
                    else
                      return(res)                                            
                  },
                  objectiveFunctionParams = list(caseClass = caseClass, saveModel = saveModel),
                  direction = "maximize",
                  name = "CV.NPV")
}

# Predefined objective calculating the confusion of two classes
# in a cross-validation experiment
cvConfusion <- function(nfold=10, ntimes=10, leaveOneOut=FALSE, stratified = FALSE, foldList=NULL,
                        trueClass, predictedClass, saveModel=FALSE)
{
  createObjective(precalculationFunction = crossValidation,
                  precalculationParams = list(nfold = nfold, 
                                              ntimes = ntimes, 
                                              leaveOneOut = leaveOneOut,
                                              stratified = stratified,
                                              foldList = foldList),
                  objectiveFunction = function(result, trueClass, predictedClass, saveModel)
                                      {
                                        res <- mean(sapply(result,
                                                 function(run)
                                                      {
                                                        predictedLabels <- unlist(lapply(run,
                                                                                  function(fold)fold$predictedLabels))
                                                        trueLabels <- unlist(lapply(run,
                                                                             function(fold)fold$trueLabels))
                                                        
                                                        return(sum(predictedLabels[trueLabels == trueClass]
                                                                   == predictedClass,na.rm=TRUE) / 
                                                               sum(trueLabels == trueClass))
                                                      }))
                                          if (saveModel)
                                          {
                                            model <- lapply(result,function(run)
                                                      {
                                                        lapply(run,function(fold)
                                                        {
                                                          fold$model
                                                        })
                                                      })
                                            return(list(additionalData=model, fitness=res))
                                          }
                                          else
                                            return(res)
 
                                      },
                  objectiveFunctionParams = list(trueClass=trueClass, predictedClass=predictedClass, saveModel=saveModel),
                  direction="minimize",
                  name="CV.Confusion")
}
