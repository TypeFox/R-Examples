## functions to train and test random forest

trainRF = function(labelDir, featDirs, names, combineStanding=FALSE, strat=TRUE, 
                   ntree=500, mtry=NULL, replace=TRUE,
                   nsample=10000, nodesize=1, sampsize=10000) {
  # function to train a random forest
  
  cat("loading training data\n")
  train = loadData(labelDir, featDirs, names)
  trainDat = train[[2]]
  trainDat$timestamp = NULL
  trainDat$PtID = NULL
  
  # pre-process - center and scale features
  cat("pre-processing\n")
  preProcValues = preProcess(trainDat, method = c("center", "scale"))
  trainDat = predict(preProcValues, trainDat)
  trainDat[is.na(trainDat)] = 0
  
  cat("training RF model\n")
  labels = train[[1]]$behavior
  trainDat = trainDat[labels!="NULL", ]
  labels = labels[labels!="NULL"]
  
  cat("using labels: ", unique(labels), "\n")
  
  if (!is.null(nsample)){
    nsample = min(nrow(trainDat), nsample)
  } else { nsample = nrow(trainDat) }
  if (strat) {
    nstrat = round(nsample / length(unique(labels)))
    s = stratSample(labels, nstrat)
  } else {
    s = sample(nrow(trainDat), nsample)
  }
  sampsize = min(length(s), sampsize)
  if (is.null(mtry)) { mtry = floor(sqrt(ncol(trainDat))) }
  rf = randomForest(trainDat[s, ], factor(labels[s]),
                    ntree=ntree,
                    mtry=mtry,
                    replace=replace,
                    sampsize=sampsize,
                    nodesize=nodesize,
                    importance=TRUE)
  rf$preProcValues = preProcValues
  rf$groundTruth = labels[s]
  return(rf)
}
testRF = function(featDirs, modelName, saveDir, testNames){
  # load model
  rf = loadModel(modelName, "rf")
  
  if (length(testNames) == 0) {
    stop("no test data\n")
  }
  for (testName in testNames) {
    testDat = loadFeatures(featDirs, testName)
    if (nrow(testDat) == 0) {
      cat(file.path(featDirs, testName), "not found\n")
      next
    }
      
    cat(testName, "\n")
    # remove timestamps
    timestamps = testDat$timestamp
    testDat$timestamp = NULL
    
    # pre-process - center and scale features
    cat("pre-processing\n")
    cat(length(names(testDat)), "\n")
    cat(rf$preProcValues$dim[2], "\n")
    if (length(names(testDat)) != rf$preProcValues$dim[2]) {
      stop("test data dimensions don't match model dimensions")
    }
    cat("dimensions matched\n")
    testDat = predict(rf$preProcValues, testDat)
    cat("preprocessed\n")
    testDat[is.na(testDat)] = 0
    
    cat("applying RF model\n")
    pr = predict(rf, testDat)
    saveFile = file.path(saveDir, paste0(testName,".csv"))
    writePredictions(pr, timestamps, saveFile)
  }
}
testTwoRFs = function(featDirs1, featDirs2, rf1, rf2, saveDir, testNames){
  if (length(testNames) == 0) {
    stop("no test data\n")
  }
  for (testName in testNames) {
    
    # first rf
    testDat1 = loadFeatures(featDirs1, testName)
    if (nrow(testDat1) == 0) {
      cat(file.path(featDirs1, testName), "not found\n")
      next
    }
    testDat2 = loadFeatures(featDirs2, testName)
    if (nrow(testDat2) == 0) {
      cat(file.path(featDirs2, testName), "not found\n")
      next
    }
    
    cat(testName, "\n")
      
    # remove timestamps
    timestamps = testDat1$timestamp
    testDat1$timestamp = NULL
    # pre-process - center and scale features
    cat("pre-processing\n")
    if (length(names(testDat1)) != rf1$preProcValues$dim[2]) {
      stop("test data dimensions don't match model dimensions")
    }
    testDat1 = predict(rf1$preProcValues, testDat1)
    testDat1[is.na(testDat1)] = 0
    cat("applying RF model\n")
    pr1 = predict(rf1, testDat1, type="prob")
    pr1 = data.frame(pr1)
    pr1$timestamp = timestamps
    pr1 = pr1[ ,c("timestamp", "Sedentary", "StandingStill", "StandingMoving",
                    "Walking", "Biking", "Vehicle")]
  
    # second rf
    # remove timestamps
    timestamps = testDat2$timestamp
    testDat2$timestamp = NULL
    # pre-process - center and scale features
    #cat("pre-processing\n")
    if (length(names(testDat2)) != rf2$preProcValues$dim[2] ){
      stop("test data dimensions don't match model dimensions")
    }
    testDat2 = predict(rf2$preProcValues, testDat2)
    testDat2[is.na(testDat2)] = 0
    #cat("applying RF model\n")
    pr2 = predict(rf2, testDat2, type="prob")
    pr2 = data.frame(pr2)
    pr2$timestamp = timestamps
    
    pr2$Biking = pr2$Cycling
    pr2$Vehicle = pr2$CarRiding
    pr2$Sedentary = 2*pr2$Sitting/3
    pr2$Walking = pr2$Walking + pr2$Jogging
    pr2$StandingMoving = pr2$Sitting/6
    pr2$StandingStill = pr2$Sitting/6
    pr2 = pr2[ ,c("timestamp", "Sedentary", "StandingStill", "StandingMoving", 
                  "Walking", "Biking", "Vehicle")]
    
    pr1 = pr1[pr1$timestamp %in% pr2$timestamp, ]
    pr2 = pr2[pr2$timestamp %in% pr1$timestamp, ]
    timestamps = pr1$timestamp
    pr1$timestamp = NULL
    pr2$timestamp = NULL
    
    pr = pr1 + pr2
    prr = character(0)
    for (i in 1:nrow(pr)) {
      prr = c(prr, names(which.max(pr[i, ])))
    }
    saveFile = file.path(saveDir, paste0(testName, ".csv"))
    writePredictions(prr, timestamps, saveFile)
  }
}
stratSample = function(labels, nsamp) {
  t = table(labels)
  lmat = data.frame(label=labels)
  lmat$idx = as.numeric(rownames(lmat))
  samples = numeric(0)
  for (i in 1:length(names(t))) {
    s = sample(lmat[lmat$label == names(t)[i], c("idx")], nsamp, replace=TRUE)
    samples = c(samples, s)
  }
  return(samples)
}