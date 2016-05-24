## functions to load and save data & other extras

loadModel = function(modelName, which) {
  # function to load a model
  hmm = NULL
  rf = NULL
  winSize = NULL
  
  if (typeof(modelName) == "list") {
    return(modelName)
  }
  if (file.exists(modelName)) {
    # if the model is a path to file, load it
    load(modelName)
  } else {
    # otherwise look in package data
    data(list=modelName, package="classifyBehaviors", envir = environment())
  }
  if (which=="winSize") {
    return(winSize)
  } else if (which=="rf") {
    return(rf)
  } else if (which=="hmm") {
    return(hmm)
  }
}
loadData = function(labelDir, featDirs, names=NULL) {
  cat(featDirs,"\n")
  if (is.null(names)) {
    names = list.files(labelDir)
  }
  if (length(names) == 0) {
    return(NULL)
  }
  all_data = data.frame()
  all_labels = data.frame()
  for (name in names) {
    cat(name, "\n")
    # for each participant
    days = list.files(file.path(labelDir, name))
    for (day in days) {
      # for each day
      checkFlag = TRUE
      cat(" ", file_path_sans_ext(day), "\n")
      for (featDir in featDirs) {
        # for each feature type
        featFile = file.path(featDir, name, file_path_sans_ext(day))
        if (!file.exists(featFile)) {
          #skip this day
          checkFlag = FALSE
          break
        }
      }
      if (checkFlag) {
        labels = read.csv(file.path(labelDir, name, day), stringsAsFactors=FALSE)
        data = data.frame(timestamp=labels$timestamp, stringsAsFactors=TRUE)
        for (featDir in featDirs) {
          # for each feature type
          featFile = file.path(featDir, name, file_path_sans_ext(day))
          d = read.csv(featFile, header=TRUE, stringsAsFactors=FALSE)
          data = merge(d, data, by = "timestamp", all=FALSE)
        }
        l2 = data.frame(timestamp=data$timestamp, stringsAsFactors=FALSE)
        labels = merge(labels, l2, by="timestamp", all=FALSE)
        if (nrow(data) == nrow(labels)) {
          all_data = rbind(all_data, data)
          all_labels = rbind(all_labels, labels)
        } else {
          stop("error:", name, day, "\n")
        } 
      }
    }
  }
  return(list(all_labels, all_data))
}
loadFeatures = function(featDirs, names=NULL) {
  # function to load features from directories
  if (is.null(names)) {
    names = list.files(featDirs[1])
  }
  if (length(names) == 0) {
    return(NULL)
  }
  all_feats = data.frame()
  for (i in 1:length(names)) {
    days = list.files(file.path(featDirs[1], names[i]))
    if (length(days) == 0) {
      next
    }
    for (k in 1:length(days)) {
      # for each day
      checkFlag = TRUE
      if (length(featDirs) > 1) {
        for (j in 2:length(featDirs)) {
          # for each feature type
          featFile = file.path(featDirs[j], names[i], days[k])
          if (!file.exists(featFile)) {
            #skip this day
            checkFlag = FALSE
            break
          }
        }
      }
      if (checkFlag) {
        featFile = file.path(featDirs[1], names[i], days[k])
        feats = read.csv(featFile, header=TRUE, stringsAsFactors=FALSE)
        if (length(featDirs) > 1) {
          for (j in 2:length(featDirs)) {
            # for each feature type
            featFile = file.path(featDirs[j], names[i], days[k])
            f = read.csv(featFile, header=TRUE, stringsAsFactors=FALSE)
            feats = merge(f, feats, by = "timestamp", all=FALSE)
          }
        }
        all_feats = rbind(all_feats, feats)
      }
    }
  }
  return(all_feats)
}
loadPredictionsAndLabels = function(labelDir, predDir, names=NULL) {
  if (is.null(names)) {
    names = list.files(predDir)
  }
  all_predictions = data.frame()
  for (i in 1:length(names)) {
    # for each participant
    name = file_path_sans_ext(names[i])
    predFile = paste0(file.path(predDir, name), ".csv")
    labelFile = file.path(labelDir, name)
    
    if (file.exists(predFile) & (file.exists(labelFile) | file.exists(paste0(labelFile, ".csv")))){
      predictions = read.csv(predFile, stringsAsFactors=FALSE)
      labels = loadLabels(labelDir, name)
      predictions = merge(labels, predictions, by="timestamp", all=FALSE)
      all_predictions = rbind(all_predictions, predictions)
    }
  }
  return(all_predictions)
}
loadLabels = function(labelDir, names=NULL) {
  if (is.null(names)) {
    names = list.files(labelDir)
  }
  all_labels = data.frame()
  for (i in 1:length(names)) {
    # for each participant
    if (file.exists(file.path(labelDir, names[i]))){
      labelFiles = list.files(file.path(labelDir, names[i]))
      for (k in 1:length(labelFiles)) {
        # for each day
        labels = read.csv(file.path(labelDir, names[i], labelFiles[k]), 
                          stringsAsFactors=FALSE)
        labels$id = names[i]
        all_labels = rbind(all_labels, labels)
        }
    } else {
      labels = read.csv(file.path(labelDir, names[i]), stringsAsFactors=FALSE)
      labels$id = names[i]
      all_labels = rbind(all_labels, labels)
    }
  }
  return(all_labels)
}
loadPredictions = function(predDir, names=NULL) {
  if (is.null(names)) {
    names = list.files(predDir)
  }
  all_predictions = data.frame()
  for (i in 1:length(names)) {
    # for each participant
    name = file_path_sans_ext(names[i])
    predFile = file.path(predDir, paste0(name, ".csv"))
    if (file.exists(predFile)) {
      predictions = read.csv(predFile, stringsAsFactors=FALSE)
      predictions$id = name
      all_predictions = rbind(all_predictions, predictions)
    }
  }
  return(all_predictions)
}
writePredictions = function(values, timestamps, saveFile) {
  if (length(values) != length(timestamps)) {
    stop("lengths dont match")
  }
  if (file.exists(saveFile)) {
    warning("overwriting file")
    file.remove(saveFile)
  }
  if (!file.exists(dirname(saveFile))) {
    dir.create(dirname(saveFile), recursive=TRUE)
  }
  cat("timestamp,prediction\n", file=saveFile, sep="", append=TRUE)
  for (i in 1:length(values)) {
    str = paste(timestamps[i], as.character(values[i]), sep = ",")
    cat(str, "\n", file=saveFile, sep = "", append=TRUE)
  }
}
clearFiles = function(dir) {
  # function to delete files in a directory
  if (!file.exists(dir)) {
    stop("directory doesn't exist")
  }
  if (!file.info(dir)$isdir) {
    stop("isn't a directory")
  }
  files = list.files(dir)
  for (file in files) {
    file.remove(file.path(dir, file))
  }
}
getDateFmt = function(inputString) {
  dF1 = "%Y-%m-%d %H:%M:%S"
  dF2 = "%m/%d/%Y %H:%M:%S"
  if (!is.na(strptime(str_trim(inputString), dF1))) {
    return(dF1)
  }
  if (!is.na(strptime(str_trim(inputString), dF2))) {
    return(dF2)
  }
  return(NULL)
}
alignStart = function(winSize, start) {
  d0 = trunc(start, "days")
  s = as.numeric(difftime(start, d0, units="secs"))
  w = ceiling(s / winSize)
  newStart = as.POSIXlt(d0 + w * winSize)
  return(newStart)
}