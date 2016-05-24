fscaret<-function(trainDF, testDF, installReqPckg=FALSE,
		  preprocessData=FALSE, with.labels=TRUE, classPred=FALSE,
		  regPred=TRUE, skel_outfile=NULL,
		  impCalcMet="RMSE&MSE", myTimeLimit=24*60*60,
		  Used.funcRegPred=NULL, Used.funcClassPred=NULL,
		  no.cores=NULL, method="boot", returnResamp="all",
		  missData=NULL, supress.output=FALSE, saveModel=FALSE, ... ){


mySystem <- .Platform$OS.type
regPredRES <- list()
classPredRES <- list()
impCalcRES <- list()
fscaretRES <- list()
methodSet <- method
returnResampSet <- returnResamp
fitControl <- trainControl(method = methodSet, returnResamp = returnResampSet, ...)
no.cores<-no.cores


# Get working dir
mywd <- getwd()

# Set working dir to tempdir
setwd(tempdir())

# Dummy object
data(requiredPackages, envir=environment())

# prevent Java from Java requested System.exit(130), closing R. so that JVM doesn't steal SIGINT from R 
options(java.parameters="-Xrs")

options(warn=-1)

# prevent gsubfn from loading tckl which causes problem with mclapply
options(gsubfn.engine = "R")
		  
if(installReqPckg==TRUE){

# Install required packages
try(installPckg(requiredPackages))

# Try loading required packages
loadedPackages <- try(for (i in 1:length(requiredPackages)){
require(requiredPackages[i], quietly = FALSE,
        character.only = TRUE)})

} else {

# Try loading required packages
loadedPackages <- try(for (i in 1:length(requiredPackages)){
require(requiredPackages[i], quietly = FALSE,
        character.only = TRUE)})
        
        }

# Inform user if every package was loaded successfully        
if(class(loadedPackages) != "try-error"){
cat("\n----Packages loaded successfully----\n")
cat("\n")

} else {

cat("\n----Loading required packages failed----\n")
cat("\n----Please check if you have installed:----\n")
cat("\n",requiredPackages, "\n")
stop()

}

# Check if loaded data is dataFrame
if(!(is.data.frame(trainDF))){
cat("\n----Provided data is not data.frame object----\n")
cat("\n----Please check the result of: is.data.frame(yourData) function ----\n")
}

if(!(is.data.frame(testDF))){
cat("\n----Provided data is not data.frame object----\n")
cat("\n----Please check the result of: is.data.frame(yourData) function ----\n")
}



# Set local settings back to "normal", because loading RWeka changes locale settings
Sys.setlocale(category = "LC_NUMERIC", locale = "C")

# Set models data set to use in funcRegPred
if(regPred==TRUE){

if(is.null(Used.funcRegPred)){

definedModels <- c("rf")

} else if(Used.funcRegPred=="all"){

definedModels <- funcRegPred

} else {

definedModels <- Used.funcRegPred

}

} else if(is.null(regPred)){
  
  regPred <- FALSE
  
} 

# Set models data set to use in funcClassPred
if(classPred==TRUE){

if(is.null(Used.funcClassPred)){

definedModels <- c("rf")

} else if(Used.funcClassPred=="all"){

definedModels <- funcClassPred

} else {

definedModels <- Used.funcClassPred

}

} else if(is.null(classPred)){
  
  classPred <- FALSE
  
} 

# Check for NULL skel_outfile obj

if(is.null(skel_outfile)){
  
 skel_outfile <- c("_default_") 
  
}

# Check the number of selected cores - if NULL use all available or set no.cores=1 on Windows

if(is.null(no.cores)){
  
  no.cores<-detectCores()
  
} else {
  
  no.cores <- no.cores
  
}


# if(mySystem!="windows"){
#   
#   if(is.null(no.cores)){
#   
#     no.cores<-detectCores()
#   
#   } else {
#   
#   no.cores <- no.cores
#   
#   }
#   
# } else {
# 
# no.cores <- 1
# 
# }


if(regPred==TRUE){
# Scan dimensions of trainDF [lk_row x lk_col]
lk_col = ncol(trainDF)
lk_row = nrow(trainDF)

# Read labels of trainDF
labelsFrame <- as.data.frame(colnames(trainDF[1:(ncol(trainDF)-1)]))

# Create a train data set matrix
trainMatryca_nr <- matrix(data=NA,nrow=lk_row,ncol=lk_col)
colnames(trainMatryca_nr) <- colnames(trainDF)

row=0
col=0

for(col in 1:(lk_col)) {
   for(row in 1:(lk_row)) {
     trainMatryca_nr[row,col] <- (as.double(trainDF[row,col]))
    }
}


# Scan dimensions of trainDF [lk_row x lk_col]
lk_col_test = ncol(testDF)
lk_row_test = nrow(testDF)

testMatryca_nr <- matrix(data=NA,nrow=lk_row_test,ncol=lk_col_test)
colnames(testMatryca_nr) <- colnames(testDF)

row=0
col=0

for(col in 1:(lk_col_test)) {
   for(row in 1:(lk_row_test)) {
     testMatryca_nr[row,col] <- (as.double(testDF[row,col]))
    }
  }

}

if(classPred==TRUE){

# Scan dimensions of trainDF [lk_row x lk_col]
lk_col = ncol(trainDF)
lk_row = nrow(trainDF)

# Read labels of trainDF
labelsFrame <- as.data.frame(colnames(trainDF[1:(ncol(trainDF)-1)]))

# Create a train data set matrix
trainMatryca_nr <- trainDF
colnames(trainMatryca_nr) <- colnames(trainDF)

# Scan dimensions of trainDF [lk_row x lk_col]
lk_col_test = ncol(testDF)
lk_row_test = nrow(testDF)

testMatryca_nr <- testDF
colnames(testMatryca_nr) <- colnames(testDF)

}



# Check for missing data
if(!is.null(missData)){
  if(missData=="delRow"){

# record rows with missing values
missing.rowsTrain <- which(rowSums(is.na(trainMatryca_nr))>0)
missing.rowsTest <- which(rowSums(is.na(testMatryca_nr))>0)

if(length(missing.rowsTrain)>0){

# Delete rows with at least one missing value
trainMatryca_nr <- trainMatryca_nr[-missing.rowsTrain,]
lk_row <- nrow(trainMatryca_nr)

# Show warning message
cat("\n", "Warning!","\n" ,"Rows:","\n")
print(as.data.frame(missing.rowsTrain))
cat("\n", " from training data set were deleted because of missing values.","\n")

}

if(length(missing.rowsTest)>0){

# Delete rows with at least one missing value
testMatryca_nr <- testMatryca_nr[-missing.rowsTest,]
lk_row_test <- nrow(testMatryca_nr)

# Show warning message
cat("\n","Warning!","\n" ,"Rows:","\n")
print(as.data.frame(missing.rowsTest))
cat("\n", " from testing data set were deleted because of missing values.","\n")

    }
  } 

  if(missData=="delCol"){

tmpMatrix <- rbind(trainMatryca_nr, testMatryca_nr)

# record cols with missing values
missing.colsTmpMatrix <- which(colSums(is.na(tmpMatrix))>0)

if(length(missing.colsTmpMatrix)>0){

# Delete cols with at least one missing value
trainMatryca_nr <- trainMatryca_nr[,-missing.colsTmpMatrix]
testMatryca_nr <- testMatryca_nr[,-missing.colsTmpMatrix]
labelsFrame <- subset(labelsFrame, select=-missing.colsTmpMatrix)

lk_col <- ncol(trainMatryca_nr)
lk_col_test <- ncol(testMatryca_nr)

# Show warning message
cat("\n","Warning!","\n" ,"Cols:","\n")
print(as.data.frame(missing.colsTmpMatrix))
cat("\n", " from training and testing data set were deleted because of missing values.","\n")

    }
  }

  if(missData=="meanCol"){

# Bind matricies
tmpMatrix <- rbind(trainMatryca_nr, testMatryca_nr)

# record cols with missing values
missing.colsTmpMatrix <- which(colSums(is.na(tmpMatrix))>0)

if(length(missing.colsTmpMatrix)>0){

# Show warning message
cat("\n","Warning!","\n" ,"There were cols with missing data:","\n")
print(as.data.frame(missing.colsTmpMatrix))
cat("\n","Replacing NA's with column mean","\n")

# Replace missing values with column median
for(i in 1:length(missing.colsTmpMatrix)){

rowToReplace <- missing.colsTmpMatrix[i]
tmpMatrix[,rowToReplace] <- impute.mean(tmpMatrix[,rowToReplace])

}

# Unbind matricies
trainMatryca_nr <- tmpMatrix[1:(lk_row),]
testMatryca_nr <- tmpMatrix[(lk_row+1):(nrow(tmpMatrix)),]

    }
  }
}


# Data preprocess

if(preprocessData==TRUE){

preprocessRes <- dataPreprocess(trainMatryca_nr, testMatryca_nr, labelsFrame, lk_col, lk_row, with.labels)

lk_col = ncol(preprocessRes$trainMatryca)
lk_row = nrow(preprocessRes$trainMatryca)
lk_col_test = ncol(preprocessRes$testMatryca)
lk_row_test = nrow(preprocessRes$testMatryca)

trainDF <- preprocessRes$trainMatryca
testDF <- preprocessRes$testMatryca

trainMatryca_nr <- matrix(data=NA,nrow=lk_row,ncol=lk_col)
colnames(trainMatryca_nr) <- colnames(trainDF)

row=0
col=0

for(col in 1:(lk_col)) {
   for(row in 1:(lk_row)) {
     trainMatryca_nr[row,col] <- (as.double(trainDF[row,col]))
    }
}


testMatryca_nr <- matrix(data=NA,nrow=lk_row_test,ncol=lk_col_test)
colnames(testMatryca_nr) <- colnames(testDF)

row=0
col=0

for(col in 1:(lk_col_test)) {
   for(row in 1:(lk_row_test)) {
     testMatryca_nr[row,col] <- (as.double(testDF[row,col]))
    }
}

labelsFrame <- as.data.frame(preprocessRes$labelsDF)

} else if(preprocessData==FALSE){
  
  orig_input_no <- data.frame("Orig_labels"=c(1:(lk_col-1)))
  exportlabelsFrame <- data.frame(orig_input_no,labelsFrame[,1])
  colnames(exportlabelsFrame)<-c("Orig Input No","Labels")
  
  labelsFrame <- as.data.frame(exportlabelsFrame)
  
  }


# Engine for regression
if(regPred==TRUE){

# Suppress warnings
options(warn=-1)

cat("\n-----Warnings have been supressed!----\n")

# definition of input and output vector
xTrain <- data.frame(trainMatryca_nr[,-lk_col])
yTrain <- as.vector(trainMatryca_nr[,lk_col])

xTest <- data.frame(testMatryca_nr[,-lk_col])
yTest <- as.vector(testMatryca_nr[,lk_col])

regPredRES <- regVarImp(definedModels, xTrain, yTrain, xTest, fitControl,
			myTimeLimit, no.cores, lk_col, supress.output, mySystem)

			
if(is.null(impCalcMet)){

print("You haven't chosen impCalcMet, so no variable importance calculations were done!")

} else if((!is.null(impCalcMet))&&((impCalcMet=="RMSE")||(impCalcMet=="MSE")||(impCalcMet=="RMSE&MSE"))){

impCalcRES <- impCalc(skel_outfile, xTest, yTest, lk_col,labelsFrame, with.labels, regPred, classPred, saveModel)

}

fscaretRES <- list(ModelPred=regPredRES, VarImp=impCalcRES)

if(preprocessData==TRUE){
  
  fscaretRES <- list(ModelPred=regPredRES, VarImp=impCalcRES,
		      PPlabels=labelsFrame, PPTrainDF=trainDF,
		      PPTestDF=testDF)
  
}			

}

# Engine for classification
if(classPred==TRUE){

# Suppress warnings
options(warn=-1)

cat("\n-----Warnings have been supressed!----\n")

# # definition of input and output vector
xTrain <- data.frame(trainMatryca_nr[,-lk_col])
yTrain <- as.factor(trainMatryca_nr[,lk_col])

xTest <- data.frame(testMatryca_nr[,-lk_col])
yTest <- as.factor(testMatryca_nr[,lk_col])

classPredRES <- classVarImp(definedModels, xTrain, yTrain, xTest, fitControl,
			myTimeLimit, no.cores, lk_col, supress.output, mySystem)

			
if(is.null(impCalcMet)){

print("You haven't chosen impCalcMet, so no variable importance calculations were done!")

} else if((!is.null(impCalcMet))&&((impCalcMet=="RMSE")||(impCalcMet=="MSE")||(impCalcMet=="RMSE&MSE"))){

impCalcRES <- try(impCalc(skel_outfile, xTest, yTest, lk_col,labelsFrame, with.labels, regPred, classPred, saveModel))

if(class(impCalcRES)=="try-error"){

print("Importance scaling has failed. Please check Rout file. For more information set 'supress.output=FALSE'")

}

}

fscaretRES <- list(ModelPred=classPredRES, VarImp=impCalcRES)

if(preprocessData==TRUE){
  
  fscaretRES <- list(ModelPred=classPredRES, VarImp=impCalcRES,
		      PPlabels=labelsFrame, PPTrainDF=trainDF,
		      PPTestDF=testDF)
  
}			

}
# Return to your working dir
setwd(mywd)

return(fscaretRES)

}