#File that defines some KEEL Util functions

#Run a set of algorithms in parallel
#Arguments could be a list with the algorithms or various algorithms as various arguments
runParallel <- function(algorithmList, cores) {

  if(missing(cores)) {
    #Assign num of cores
    cores <- parallel::detectCores()
  }

  if(class(algorithmList) != "list") {
    stop("Error. Argument must be a list with the algorithm objects.")
  }

  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  cat(paste0("Executing experiment in " ,  foreach::getDoParWorkers(), " cores"), sep="\n")

  #Execute algorithms in parallel
  i<-1
  exec_algs <- foreach(i=1:length(algorithmList)) %dopar% {
    algorithmList[[i]]$run()
    return(algorithmList[[i]])
  }
  #Stop cluster
  parallel::stopCluster(cl)

  cat("Experiment executed successfully", sep="\n")

  return(exec_algs)
}

#Execute a set of algorithms sequentially
#Arguments could be a list with the algorithms or various algorithms as various arguments
runSequential <- function(algorithmList) {

  if(class(algorithmList) != "list") {
    stop("Error. Argument must be a list with the algorithm objects.")
  }

  cat("Executing experiment sequential", sep="\n")

  #Eexcute algorithms sequentially
  i<-1
  exec_algs <- foreach(i=1:length(algorithmList)) %do% {
    algorithmList[[i]]$run()
    return(algorithmList[[i]])
  }

  cat("Experiment executed successfully", sep="\n")

  return(exec_algs)
}


#Read dataset in keel format
#Returns a data.frame with the data
read.keel <- function(file){

  text <- readLines(file)

  i <- 1
  while(!grepl("@attribute", tolower(text[i]))){
    i <- i+1
  }
  #now text is first @attribute
  attributeNames <- c()
  attributeTypes <- c()

  while(grepl("@attribute", tolower(text[i]))){
    #Obtain attribute i name
    attributeNames <- c(attributeNames, gsub("\'", "", strsplit(text[i], "[ {]")[[1]][2]))

    if(grepl("\\{", text[i])){
      #If line contains "{", attribute is categorical
      attributeTypes <- c(attributeTypes, "categorical")
    }
    else{
      #real or integer attribute
      attributeTypes <- c(attributeTypes, strsplit(text[i], " ")[[1]][3])
    }

    i <- i+1
  }



  outputs <- -1

  while(!grepl("@data", tolower(text[i]))){
    if(grepl("@outputs", tolower(text[i]))){
      outputAttribute <- gdata::trim(strsplit(text[i], " ")[[1]][2])

      for(j in 1:length(attributeNames)){
        if(grepl(outputAttribute, attributeNames[j])){
          outputs <- j
        }
      }

      if(outputs == -1){
        stop("Output attribute don't found")
      }

    }
    i <- i+1
  }
  i <- i+1
  #now text is first data line

  data <- c()

  #num of rows
  row <- 0

  while(!is.na(text[i])){
    #Split data by commas
    dataWords <- strsplit(text[i], ",")[[1]]

    #If words are <= 0, there are no more data
    if(length(dataWords) > 0){
      dataLine <- c()

      #Create data line depending on attribute type
      for(j in 1:length(dataWords)){
        if(dataWords[j] == '?'){
          #stop("NA")
          dataLine <- c(dataLine, "NA")
        }
        else if(grepl("integer", attributeTypes[j])){
          dataLine <- c(dataLine, strtoi(gdata::trim(dataWords[j])))
        }
        else if(grepl("real", attributeTypes[j])){
          dataLine <- c(dataLine, as.double(gdata::trim(dataWords[j])))
        }
        else if(grepl("categorical", attributeTypes[j])){
          dataLine <- c(dataLine, gdata::trim(dataWords[j]))
        }
        else{
          stop("Type not found")
        }
      }

      #Add data line to full data
      data <- c(data, dataLine)
      row <- row+1
    }

    i <- i+1

  }

  #Create data matrix
  m <- matrix(data, nrow = row, ncol=length(attributeNames), byrow = TRUE)
  #Set column names
  colnames(m) <- attributeNames

  #If output is not last attribute, change columns (in matrix and attribute names and types)
  if((outputs != -1) && (outputs < length(attributeNames))){
    m2 <- m
    attributeTypes2 <- attributeTypes
    attributeNames2 <- attributeNames

    j <- 1
    for(i in 1:length(attributeNames)){
      if(i != outputs){
        m[,j] <- m2[,i]
        attributeTypes[j] <- attributeTypes2[i]
        attributeNames[j] <- attributeNames2[i]
        j <- j+1
      }
    }
    m[,length(attributeNames)] <- m2[,outputs]
    attributeTypes[length(attributeNames)] <- attributeTypes2[outputs]
    attributeNames[length(attributeNames)] <- attributeNames2[outputs]
  }

  #Generate data.frame
  df <- data.frame(m)

  #Convert categorical data to character
  for(i in 1:length(attributeTypes)){
    if(attributeTypes[i] == "categorical"){
      df[,i] <- as.character(df[,i])
    }
  }

  #Class column as factor
  df[ncol(df)] <- as.factor(df[ncol(df)][[1]])

  return(df)
}


#Load a dataset from keel repository
#Returns a data.frame with the data
loadKeelDataset <- function(dataName){

  dataList <- RKEELdata::getKeelDatasetList()

  #At the moment, only few datasets from keel repository are attached to RKEEL
  if(dataName %in% dataList) {
    dataPath <- RKEELdata::getDataPath()

    if(substr(dataPath, nchar(dataPath), nchar(dataPath)) != "/"){
      dataPath <- paste0(dataPath, "/")
    }

    cat(dataPath)

    df <- read.keel(paste0(dataPath, dataName, ".dat"))

    return(df)
  }
  else {
    cat(paste0("Dataset ", dataName, " is not available."), sep = "\n")
    cat("Please select one of the following: ", sep="\n")
    for(d in dataList){
      cat(paste0("   ", d), sep = "\n")
    }

    return(NULL)
  }

}


#Write .dat dataset file from a data frame
writeDatFromDataframe = function(data, fileName){

  dataName <- fileName

  #Check if data is a data.frame
  if(!is.data.frame(data)){
    stop(paste0("Error. Must give a data.frame."))
  }

  #full dataset string
  text <- ""

  #add relationName
  text <- paste0(text, "@relation ", dataName, "\n")

  attributesType <- c()
  #add attributes name and type
  for(i in 1:length(colnames(data))){

    #add "@attribute" and attribute name
    attribute <- paste0("@attribute ", colnames(data)[i])

    #caterogical
    if((typeof(data[,i]) == "character") || ( !is.na(match(TRUE, is.na(suppressWarnings(as.numeric(as.character(data[,i])))))) )  ){
      #add "{" and first value
      attribute <- paste0(attribute, " {", unique(data[,i])[1])
      #Start in 2 for no comma problems; add all other values
      for(l in 2:length(unique(data[,i]))){
        attribute <- paste0(attribute, ", ", unique(data[,i])[l])
      }
      #finish with "}"
      attribute <- paste0(attribute, "}")
      attributesType <- c(attributesType, "character")
    }
    #real
    else if(typeof(as.numeric(as.character(data[,i]))) == "double"){
      #add type, min and max values
      minValue <- format(min(na.omit(as.numeric(as.character(data[,i])))), nsmall = 1)
      maxValue <- format(max(na.omit(as.numeric(as.character(data[,i])))), nsmall = 1)
      attribute <- paste0(attribute, " real [", minValue, ", ", maxValue, "]")
      attributesType <- c(attributesType, "real")
    }
    #integer
    else if(typeof(as.numeric(as.character(data[,i]))) == "integer"){
      #add type, min and max values
      attribute <- paste0(attribute, " integer [", min(na.omit(as.numeric(as.character(data[,i])))), ", ", max(na.omit(as.numeric(as.character(data[,i])))), "]")
      attributesType <- c(attributesType, "integer")
    }
    #Categorical
    else if(!is.null(levels(data[,i]))){
      #add "{" and first value
      attribute <- paste0(attribute, " {", levels(data[,i])[1])
      #Start in 2 for no comma problems; add all other values
      for(l in 2:length(levels(data[,i]))){
        attribute <- paste0(attribute, ", ", levels(data[,i])[l])
      }
      #finish with "}"
      attribute <- paste0(attribute, "}")
      attributesType <- c(attributesType, "character")
    }

    #Add attribute line to full dataset string
    text <- paste0(text, attribute, "\n")
  }

  #Add "@data"
  text <- paste0(text, "@data", "\n")

  #Add data lines
  for(i in 1:nrow(data)){

    dataLine <- ""

    for(j in 1:ncol(data)){
      #add values separated with commas
      if(is.na(data[i,j]) || is.nan(data[i,j]) || is.null(data[i,j])) {
        dataLine <- paste0(dataLine, "<null>, ")
      }
      else{
        if(attributesType[j] == "real"){
          dataLine <- paste0(dataLine, format(as.numeric(as.character(data[i,j])), nsmall = 1), ", ")
          #dataLine <- paste0(dataLine, data[i,j], ", ")
        }
        else{
          dataLine <- paste0(dataLine, data[i,j], ", ")
        }

      }
    }

    #Delete last comma
    dataLine <- gsub(", $", "", dataLine)
    #Add data line to full dataset string
    text <- paste0(text, dataLine, "\n")
  }

  #Save dataset file
  fileConn<-file(fileName)
  writeLines(text, fileConn)
  close(fileConn)

}


#Write train and test .dat dataset
#Writing both together, there are fewer problems in min and max limits, and in classes
writeDatFromDataframes = function(trainData, testData, trainFileName, testFileName){

  #Check if data is a data.frame
  if((!is.data.frame(trainData)) || (!is.data.frame(testData))){
    stop(paste0("Error. Must give a data.frame."))
  }

  #full dataset string
  textTrain <- ""
  textTest <- ""

  #add relationName
  #text <- paste0(text, "@relation ", dataName, "\n")

  list_return <- getAttributeLinesFromDataframes(trainData, testData)

  textTrain <- paste0(textTrain, "@relation ", "train", "\n")
  textTest <- paste0(textTest, "@relation ", "test", "\n")
  textTrain <- paste0(textTrain, list_return[[1]])
  textTest <- paste0(textTest, list_return[[1]])

  attributesType <- list_return[[2]]

  #TrainData
  #Add "@data"
  textTrain <- paste0(textTrain, "@data", "\n")

  #Add data lines
  for(i in 1:nrow(trainData)){

    dataLine <- ""

    for(j in 1:ncol(trainData)){
      #add values separated with commas
      if(is.na(trainData[i,j]) || is.nan(trainData[i,j]) || is.null(trainData[i,j])) {
        dataLine <- paste0(dataLine, "<null>, ")
      }
      else{
        if(attributesType[j] == "real"){
          dataLine <- paste0(dataLine, format(as.numeric(as.character(trainData[i,j])), nsmall = 1), ", ")
          #dataLine <- paste0(dataLine, data[i,j], ", ")
        }
        else{
          dataLine <- paste0(dataLine, trainData[i,j], ", ")
        }

      }
    }

    #Delete last comma
    dataLine <- gsub(", $", "", dataLine)
    #Add data line to full dataset string
    textTrain <- paste0(textTrain, dataLine, "\n")
  }

  #Save dataset file
  fileConnTrain<-file(trainFileName)
  writeLines(textTrain, fileConnTrain)
  close(fileConnTrain)


  #TestData
  #Add "@data"
  textTest <- paste0(textTest, "@data", "\n")

  #Add data lines
  for(i in 1:nrow(testData)){

    dataLine <- ""

    for(j in 1:ncol(testData)){
      #add values separated with commas
      if(is.na(testData[i,j]) || is.nan(testData[i,j]) || is.null(testData[i,j])) {
        dataLine <- paste0(dataLine, "<null>, ")
      }
      else{
        if(attributesType[j] == "real"){
          dataLine <- paste0(dataLine, format(as.numeric(as.character(testData[i,j])), nsmall = 1), ", ")
          #dataLine <- paste0(dataLine, data[i,j], ", ")
        }
        else{
          dataLine <- paste0(dataLine, testData[i,j], ", ")
        }

      }
    }

    #Delete last comma
    dataLine <- gsub(", $", "", dataLine)
    #Add data line to full dataset string
    textTest <- paste0(textTest, dataLine, "\n")
  }

  #Save dataset file
  fileConnTest<-file(testFileName)
  writeLines(textTest, fileConnTest)
  close(fileConnTest)

}


#Get attribute lines from train and test dataframes
getAttributeLinesFromDataframes = function(trainData, testData){

  data <- rbind(trainData, testData)

  text <- ""

  attributesType <- c()
  #add attributes name and type
  for(i in 1:ncol(data)){

    #add "@attribute" and attribute name
    attribute <- paste0("@attribute ", colnames(data)[i])

    #caterogical
    if((typeof(data[,i]) == "character") || ( !is.na(match(TRUE, is.na(suppressWarnings(as.numeric(as.character(data[,i])))))) )  ){
      #add "{" and first value
      attribute <- paste0(attribute, " {", unique(data[,i])[1])
      #Start in 2 for no comma problems; add all other values
      for(l in 2:length(unique(data[,i]))){
        attribute <- paste0(attribute, ", ", unique(data[,i])[l])
      }
      #finish with "}"
      attribute <- paste0(attribute, "}")
      attributesType <- c(attributesType, "character")
    }
    #real
    else if(typeof(as.numeric(as.character(data[,i]))) == "double"){
      #add type, min and max values
      minValue <- format(min(na.omit(as.numeric(as.character(data[,i])))), nsmall = 1)
      maxValue <- format(max(na.omit(as.numeric(as.character(data[,i])))), nsmall = 1)
      attribute <- paste0(attribute, " real [", minValue, ", ", maxValue, "]")
      attributesType <- c(attributesType, "real")
    }
    #integer
    else if(typeof(as.numeric(as.character(data[,i]))) == "integer"){
      #add type, min and max values
      attribute <- paste0(attribute, " integer [", min(na.omit(as.numeric(as.character(data[,i])))), ", ", max(na.omit(as.numeric(as.character(data[,i])))), "]")
      attributesType <- c(attributesType, "integer")
    }
    #Categorical
    else if(!is.null(levels(data[,i]))){
      #add "{" and first value
      attribute <- paste0(attribute, " {", levels(data[,i])[1])
      #Start in 2 for no comma problems; add all other values
      for(l in 2:length(levels(data[,i]))){
        attribute <- paste0(attribute, ", ", levels(data[,i])[l])
      }
      #finish with "}"
      attribute <- paste0(attribute, "}")
      attributesType <- c(attributesType, "character")
    }

    #Add attribute line to full dataset string
    text <- paste0(text, attribute, "\n")

  }

  return(list(text, attributesType))
}


#Check if a dataset has continuous data
hasContinuousData = function(data){

  for (column in 1:ncol(data)) {
    if(typeof(data[,1]) != "character"){
      return(TRUE)
    }
  }

  return(FALSE)
}


#Check if a dataset has multiple classes
isMultiClass = function(data){

  numClasses <- length(unique(data[,ncol(data)]))

  if(numClasses > 2){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}


#Check if a dataset has missing values
hasMissingValues = function(data){

  complete_rows <- complete.cases(data)

  if(length(unique(complete_rows)) > 1){
    #If has more than one unique values, it has true and false values for complete cases
      #So, has missing values
    return(TRUE)
  }
  else if(unique(complete_rows)[1] != TRUE){
    #If all instances has missing values, all are FALSE for complete.cases
    return(TRUE)
  }
  else{
    #Has no missing values
    return(FALSE)
  }

}
