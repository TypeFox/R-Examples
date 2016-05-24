## wrapper functions

trainModel = function(annotations, accelerometers=NULL, GPS=NULL, winSize=60, 
                      modelName, names=NULL, strat=TRUE,
                      ntree=500, mtry=NULL, replace=TRUE, 
                      nsample=10000, nodesize=1, sampsize=10000) {
  # function to train a model - either from raw data or from pre-computed features

  # annotations
  labelDir = annotationsToLabels(annotations, winSize, names)

  # features
  featDirs = sensorsToFeatures(accelerometers, GPS, winSize, names)
  if (length(featDirs) == 0) { stop("no data directories found") }
    
  #train
  cat("\ntraining model from", length(featDirs), "devices\n")
  trainFromFeatures(labelDir, featDirs, modelName=modelName, winSize=winSize, 
              names=names, strat=strat, ntree=ntree, mtry=mtry, replace=replace,
              nsample=nsample, nodesize=nodesize, sampsize=sampsize)
}
classify = function(accelerometers=NULL, GPS=NULL, modelName, saveDir, names=NULL) {
  # function to classify data - either from raw data or from pre-computed features
  
  winSize = loadModel(modelName, "winSize")
  
  # features
  featDirs = sensorsToFeatures(accelerometers, GPS, winSize, names)
  if (length(featDirs) == 0) { stop("No data directories found") }
  
  cat("\n")
  # test
  testAllDir(featDirs, modelName, saveDir, names)
}
looXval = function(annotations, accelerometers=NULL, GPS=NULL, winSize=60, 
                   saveDir, names=NULL, strat=TRUE) {
  # annotations
  labelDir = annotationsToLabels(annotations, winSize, names)
  # features
  featDirs = sensorsToFeatures(accelerometers, GPS, winSize, names)
  if (length(featDirs) == 0) { stop("no data directories found") }
    
  #train
  cat("cross-validating model from", length(featDirs), "devices\n")
  looXvalFromFeats(labelDir, featDirs, saveDir, names, strat)
}
sensorsToFeatures = function(accelerometers=NULL, GPS=NULL, winSize, names=NULL) {
  # get feature directories from raw sensor directories
  
  featDirs = character(0)
  
  # GPS
  if (!is.null(GPS)) {
    if (!file.exists(GPS)) {
      stop("GPS file/directory not found")
    }
    if (file.info(GPS)$isdir) {
      if (isFeatureDirectory(GPS)) {
        GPSFeatDir = GPS
      } else {
        GPSFeatDir = paste(GPS, "Features", as.character(winSize), sep="_")
        extractFeatsPALMSDir(GPS, GPSFeatDir, winSize, names)
      }
    } else {
      GPSFeatDir = paste(file_path_sans_ext(GPS), "Features", as.character(winSize), 
                         sep="_")
      extractFeatsPALMSOneFile(GPS, GPSFeatDir, winSize)
    }
    featDirs = c(featDirs, GPSFeatDir)
  }
  
  # accelerometers
  if (!is.null(accelerometers)) {
    for (acc in accelerometers) {
      if (!file.exists(acc)) {
        stop("accelerometer directory not found")
      }
      if (isFeatureDirectory(acc)) {
        accFeatDir = acc
      } else {
        accFeatDir = paste(acc, "Features", as.character(winSize), sep="_")
        extractAccelerometerFeatures(acc, accFeatDir, winSize, names)
      }
      featDirs = c(featDirs, accFeatDir)
    }
  }
  return(featDirs)
}
isFeatureDirectory = function(dir) {
  # check if the directory is a feature directory
  file = list.files(dir, full.names=TRUE)[1]
  if (file.info(file)$isdir) {
    return (TRUE)
  } else {
    return (FALSE)
  }
}
looXvalFromFeats = function(labelDir, featDirs, saveDir, names=NULL, strat=TRUE) {
  saveDir1 = paste(saveDir, "Temp", sep="")
  if (is.null(names)) {
    names = list.files(labelDir)
  }
  for (i in 1:length(names)) {
    cat("test subject:", names[i], "\n")
    testNames = names[i]
    trainNames = names[-i]

    # do two-level classification
    modelName1 = "temp.rf"
    modelName2 = "temp.hmm"
    
    # first train RF
    rf = trainRF(labelDir, featDirs, trainNames, strat=strat)
    testRF(featDirs, rf, saveDir1, testNames)
    # calculate performance
    cat(testNames, "\n")
    calcPerformance(featDirs, saveDir1, testNames)
    
    # then apply HMM smoothing to RF outputs
    hmm = trainHMM(labelDir, rf, trainNames)
    testHMM(saveDir1, hmm, saveDir, testNames)
    # calculate performance
    cat(testNames,"\n")
    calcPerformance(labelDir, saveDir, testNames)
    file.remove(modelName1)
    file.remove(modelName2)
    cat("----------------------------------\n")
  }
  cat("Overall RF\n")
  calcPerformance(labelDir, saveDir1, names)
  cat("Overall HMM\n")
  calcPerformance(labelDir, saveDir, names)
  cat("----------------------------------\n")
}
trainFromFeatures = function(labelDir, featDirs, modelName, winSize, names=NULL, 
                       strat=TRUE, ntree=500, mtry=NULL, sampsize=10000,
                       replace=TRUE, nsample=10000, nodesize=1) {
  if (is.null(names)) {
    names = list.files(labelDir)
  }
  rf = trainRF(labelDir, featDirs, names=names, strat=strat, ntree=ntree, 
               mtry=mtry, replace=replace, nsample=nsample, nodesize=nodesize,
               sampsize=sampsize)
  hmm = trainHMM(labelDir, rf, names)
  if (!file.exists(dirname(modelName))){
    dir.create(dirname(modelName), recursive=TRUE)
  }
  save(rf, hmm, winSize, file=modelName)
  cat("model saved to", modelName, "\n")
}
testAllDir = function(featDirs, modelName, saveDir, names=NULL) {
  if (is.null(names)) {
    names = list.files(featDirs[1])
  }
  saveDir1 = file.path(saveDir, "Temp")
  for (i in 1:length(names)) {
    testRF(featDirs, modelName, saveDir1, names[i])
    testHMM(saveDir1, modelName, saveDir, names[i])
  }
  cat("predictions saved to", saveDir, "\n")
}