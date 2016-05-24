impCalc <- function(skel_outfile, xTest, yTest, lk_col, labelsFrame, with.labels, regPred, classPred,saveModel){
# Get RMSE from all .RData files

filesRData <- list.files(pattern = "_default_.*.RData")
impCalcScaleRMSE <- list()
impCalcScaleMSE <- list()
impCalcScaleF <- list()
rawMeasure <- list()
impCalcRes <- list()
modelRes <- list()
labelsDF <- as.data.frame(labelsFrame[1:(ncol(xTest)),])



# Dummy
res <- NULL
rm(res)


if((length(filesRData) > 0)&&(regPred==TRUE)){

filesRDataCols <- gsubfn(".*Control_","",filesRData)
filesRDataCols <- gsubfn(".RData","",filesRDataCols)

varImpRDataOutfileRMSE <- paste("./",skel_outfile, "VarImpRes_RMSE.RData",sep="")
varImpTxtOutfileRMSE <- paste("./",skel_outfile, "VarImpRes_RMSE.txt",sep="")

varImpRDataOutfileMSE <- paste("./",skel_outfile, "VarImpRes_MSE.RData",sep="")
varImpTxtOutfileMSE <- paste("./",skel_outfile, "VarImpRes_MSE.txt",sep="")

cat("\n----Processing files:----\n")
print(filesRData)

for (i in 1:length(filesRData)){

load(filesRData[i])

print("")
print("Calculating error for model:")
print(filesRData[i])
print("")

nameTmp <- filesRDataCols[i]

# Save model into R obj
if(saveModel==TRUE){
      for(k in 1:length(res)){
        lstName <- names(res[k])
        modelRes[[res$method]] <- c(modelRes[[res$method]],res[k])
      
      }
}

predTmp <- try(predict(res, xTest),silent=TRUE)

  if(class(predTmp)!="try-error"){

    if(is.factor(predTmp)==TRUE){

      predTmp <- as.character(predTmp)
      predTmp <- gsubfn(",",".",predTmp)
      Sys.setlocale(category = "LC_NUMERIC", locale = "C")
      predTmp <- as.numeric(predTmp)

    }

      impCalcScaleRMSE[[i]]<-RMSE(predTmp, yTest, length(yTest))
      impCalcScaleMSE[[i]]<-MSE(predTmp, yTest, length(yTest))

  } else {

      impCalcScaleRMSE[[i]] <- NA
      impCalcScaleMSE[[i]] <- NA

  }

names(impCalcScaleRMSE)[[i]] <- nameTmp
names(impCalcScaleMSE)[[i]] <- nameTmp

}


# Set local settings back to "normal", because loading RWeka changes locale settings
Sys.setlocale(category = "LC_NUMERIC", locale = "C")

impCalcScaleRMSE <- as.data.frame(impCalcScaleRMSE)
impCalcScaleMSE <- as.data.frame(impCalcScaleMSE)

maxRmse <- try(max(impCalcScaleRMSE, na.rm=TRUE),silent=TRUE)

if(class(maxRmse)=="try-error"){
maxRmse <- 100000
}

maxMse <- try(max(impCalcScaleMSE, na.rm=TRUE),silent=TRUE)

if(class(maxMse)=="try-error"){
maxMse <- 100000
}

impCalcScaleRMSE[is.na(impCalcScaleRMSE)]<-100000*maxRmse
impCalcScaleMSE[is.na(impCalcScaleMSE)]<-100000*maxMse

rawRMSE <- as.data.frame(impCalcScaleRMSE)
rawMSE <- as.data.frame(impCalcScaleMSE)

minRmse <- min(impCalcScaleRMSE)
minMse <- min(impCalcScaleMSE)

impCalcScaleRMSE <- (impCalcScaleRMSE/minRmse)^(-1)
impCalcScaleMSE <- (impCalcScaleMSE/minMse)^(-1)


# Get all the results from txt files and put them into one dataframe

# Load list of files with _VarImp_

filesVarImp <- list.files(pattern = "_VarImp_")

matrycaVarImp.RMSE<-matrix(data=0,nrow=(lk_col-1),ncol=(length(filesVarImp)+3))
matrycaVarImp.MSE<-matrix(data=0,nrow=(lk_col-1),ncol=(length(filesVarImp)+3))


cat("\n----Processing files:----\n")
print(filesVarImp)

# Concatenate the results
for(i in 1:length(filesVarImp)){

currentFile <- filesVarImp[i]

# read file
tempDF <- read.csv(filesVarImp[i],header=TRUE,sep="\t", strip.white = TRUE, na.strings = c("NA",""))

if(ncol(tempDF) > 1){

tempDF <- as.data.frame(rowSums(tempDF))
colnames(tempDF) <- "Overall"

}

if(is.factor(tempDF[,1])==TRUE){

tempDF[,1] <- as.character(tempDF[,1])
tempDF[,1] <- gsubfn(",",".",tempDF[,1])
Sys.setlocale(category = "LC_NUMERIC", locale = "C")
tempDF[,1] <- as.numeric(tempDF[,1])

}

# Check if there are any NA values and zero them before summing
tempDF[is.na(tempDF)]<-0
# Get absolute values of variable importance
tempDF <- abs(tempDF)


if(with.labels==FALSE){

# Delete V char from rownames
rownames(tempDF)<-gsubfn("V","",rownames(tempDF))

# Sort data.frame according to rownames
tempDF<-tempDF[order(as.numeric(rownames(tempDF))),,drop=FALSE]
# tempDF.MSE<-tempDF.MSE[order(as.numeric(rownames(tempDF.MSE))),,drop=FALSE]

} else if(with.labels==TRUE){
  
rownames(labelsDF) <- labelsDF[,2]

tempDF <- as.data.frame(tempDF[rownames(labelsDF),])

}

# Scale results from 0 to 100
tempDF.RMSE<-(tempDF/sum(tempDF[,1]))*100*impCalcScaleRMSE[,i]
tempDF.MSE<-(tempDF/sum(tempDF[,1]))*100*impCalcScaleMSE[,i]


cols=i
rows=1
for(rows in 1:nrow(tempDF.RMSE)){
matrycaVarImp.RMSE[rows,cols]<-tempDF.RMSE[rows,]
}

cols=i
rows=1
for(rows in 1:nrow(tempDF.MSE)){
matrycaVarImp.MSE[rows,cols]<-tempDF.MSE[rows,]
}

}

matrycaVarImp.RMSE[is.na(matrycaVarImp.RMSE)]<-0
matrycaVarImp.MSE[is.na(matrycaVarImp.MSE)]<-0

# Sum row-by-row
for (rows in 1:nrow(matrycaVarImp.RMSE)){
matrycaVarImp.RMSE[rows,(ncol(matrycaVarImp.RMSE)-2)]<-sum(matrycaVarImp.RMSE[rows,])
}

for (rows in 1:nrow(matrycaVarImp.MSE)){
matrycaVarImp.MSE[rows,(ncol(matrycaVarImp.MSE)-2)]<-sum(matrycaVarImp.MSE[rows,])
}

# Calculate percentages
maks.sum.rmse <- max(matrycaVarImp.RMSE[,(ncol(matrycaVarImp.RMSE)-2)])
maks.sum.mse <- max(matrycaVarImp.MSE[,(ncol(matrycaVarImp.MSE)-2)])

for (rows in 1:nrow(matrycaVarImp.RMSE)){
matrycaVarImp.RMSE[rows,(ncol(matrycaVarImp.RMSE)-1)]<-(matrycaVarImp.RMSE[rows,(ncol(matrycaVarImp.RMSE)-2)])/maks.sum.rmse*100
}

for (rows in 1:nrow(matrycaVarImp.MSE)){
matrycaVarImp.MSE[rows,(ncol(matrycaVarImp.MSE)-1)]<-(matrycaVarImp.MSE[rows,(ncol(matrycaVarImp.MSE)-2)])/maks.sum.mse*100
}

matrycaVarImp.RMSE <- as.data.frame(matrycaVarImp.RMSE)
matrycaVarImp.MSE <- as.data.frame(matrycaVarImp.MSE)

filesVarImpCols <- gsubfn(".*_VarImp_","",filesVarImp)
filesVarImpCols <- gsubfn(".txt","",filesVarImpCols)

colnames(matrycaVarImp.RMSE)<- filesVarImpCols
colnames(matrycaVarImp.MSE)<- filesVarImpCols

names(matrycaVarImp.RMSE)[length(matrycaVarImp.RMSE)-2]<-"SUM"
names(matrycaVarImp.MSE)[length(matrycaVarImp.MSE)-2]<-"SUM"

names(matrycaVarImp.RMSE)[length(matrycaVarImp.RMSE)-1]<-"SUM%"
names(matrycaVarImp.MSE)[length(matrycaVarImp.MSE)-1]<-"SUM%"

names(matrycaVarImp.RMSE)[length(matrycaVarImp.RMSE)]<-"ImpGrad"
names(matrycaVarImp.MSE)[length(matrycaVarImp.MSE)]<-"ImpGrad"


# Sort data.frame according to SUM variable importance
matrycaVarImp.RMSE<-matrycaVarImp.RMSE[order(matrycaVarImp.RMSE[,length(matrycaVarImp.RMSE)-2],decreasing=TRUE),,drop=FALSE]
matrycaVarImp.MSE<-matrycaVarImp.MSE[order(matrycaVarImp.MSE[,length(matrycaVarImp.MSE)-2],decreasing=TRUE),,drop=FALSE]


# # Importance gradient
for (rows in 2:nrow(matrycaVarImp.RMSE)){

imp1 <- matrycaVarImp.RMSE[(rows-1),(ncol(matrycaVarImp.RMSE)-2)]
imp2 <- matrycaVarImp.RMSE[(rows),(ncol(matrycaVarImp.RMSE)-2)]

matrycaVarImp.RMSE[rows,ncol(matrycaVarImp.RMSE)] <- (imp1 - imp2)/imp1*100

}
 
for (rows in 2:nrow(matrycaVarImp.MSE)){

imp1 <- matrycaVarImp.MSE[(rows-1),(ncol(matrycaVarImp.MSE)-2)]
imp2 <- matrycaVarImp.MSE[(rows),(ncol(matrycaVarImp.MSE)-2)]

matrycaVarImp.MSE[rows,ncol(matrycaVarImp.MSE)] <- (imp1 - imp2)/imp1*100

}

# Add importance variables names as last column
matrycaVarImp.MSE$Input_no <- rownames(matrycaVarImp.MSE)
matrycaVarImp.RMSE$Input_no <- rownames(matrycaVarImp.RMSE)



impCalcRes <- list("rawMSE"=rawMSE, "rawRMSE"=rawRMSE,
                   "matrixVarImp.RMSE"=matrycaVarImp.RMSE,
                   "matrixVarImp.MSE"=matrycaVarImp.MSE,
                   "model"=modelRes)

}



if((length(filesRData) > 0)&&(classPred==TRUE)){

filesRDataCols <- gsubfn(".*Control_","",filesRData)
filesRDataCols <- gsubfn(".RData","",filesRDataCols)

varImpRDataOutfileF1 <- paste("./",skel_outfile, "VarImpRes_F1.RData",sep="")
varImpTxtOutfileRMSE <- paste("./",skel_outfile, "VarImpRes_F1.txt",sep="")

cat("\n----Processing files:----\n")
print(filesRData)

for (i in 1:length(filesRData)){

load(filesRData[i])

print("")
print("Calculating error measure for model:")
print(filesRData[i])
print("")

nameTmp <- filesRDataCols[i]

# Save model into R obj
if(saveModel==TRUE){
  for(k in 1:length(res)){
    lstName <- names(res[k])
    modelRes[[res$method]] <- c(modelRes[[res$method]],res[k])
    
  }
}


predTmp <- try(as.data.frame(predict(res, xTest),silent=TRUE))
measure.table <- cbind(predTmp,yTest)

# # Check if the DF is in factor format
# if((is.factor(as.data.frame(measure.table[,1])))&&(is.factor(as.data.frame(measure.table[,2])))){
# 
# A <- 2
# 
# } else {
# 
# A <- 1
# 
# }
# 
# Convert all values to factors in measure.table

measure.table[,1] <- factor(measure.table[,1])
measure.table[,2] <- factor(measure.table[,2])

# 
# # Check how many levels are there
# 
no.lvl <- nlevels(measure.table[,2])


if(no.lvl > 2){ 	# if lvls no is > 2 then go to F1 measure

class.error.tmp <- as.data.frame(ifelse(ifelse(measure.table[,1]==measure.table[,2],1,0)))
measure.table <- cbind(measure.table,class.error.tmp)

# TODO code the F-measure calculations when levels are > 2 
# As for v0.9 the factor will be 1
# 
# 
F.measure.tmp <- 1
# 
# 
# TODO code the F-measure calculations when levels are > 2

} else if(no.lvl==2){ 	# if lvls no == 2 then go to F-measure


pred.class <- try(misclassCounts(measure.table[,1], measure.table[,2]))

F.measure.tmp <- try(pred.class$metrics$F)

}

if(class(F.measure.tmp)=="try-error"){

impCalcScaleRMSE[[i]] <- NA

} else {

impCalcScaleF[[i]] <- F.measure.tmp

}

names(impCalcScaleF)[[i]] <- nameTmp

}

# Set local settings back to "normal", because loading RWeka changes locale settings
Sys.setlocale(category = "LC_NUMERIC", locale = "C")

impCalcScaleF <- as.data.frame(impCalcScaleF)

minF <- try(min(impCalcScaleF, na.rm=TRUE),silent=TRUE)

if (class(minF)=="try-error"){
minF <- 0.0000001
}

impCalcScaleF[is.na(impCalcScaleF)]<-0.0000001*minF

rawMeasure <- as.data.frame(impCalcScaleF)


# Get all the results from txt files and put them into one dataframe

# Load list of files with _VarImp_

filesVarImp <- list.files(pattern = "_VarImp_")

matrycaVarImp.F <- matrix(data=0,nrow=(lk_col-1),ncol=(length(filesVarImp)+3))


cat("\n----Processing files:----\n")
print(filesVarImp)

# Concatenate the results
for(i in 1:length(filesVarImp)){

currentFile <- filesVarImp[i]

# read file
tempDF <- read.csv(filesVarImp[i],header=TRUE,sep="\t", strip.white = TRUE, na.strings = c("NA",""))

if(ncol(tempDF) > 1){

tempDF <- as.data.frame(rowSums(tempDF))
colnames(tempDF) <- "Overall"

}

if(is.factor(tempDF[,1])==TRUE){

tempDF[,1] <- as.character(tempDF[,1])
tempDF[,1] <- gsubfn(",",".",tempDF[,1])
Sys.setlocale(category = "LC_NUMERIC", locale = "C")
tempDF[,1] <- as.numeric(tempDF[,1])

}

# Check if there are any NA values and zero them before summing
tempDF[is.na(tempDF)]<-0
# Get absolute values of variable importance
tempDF <- abs(tempDF)

if(with.labels==FALSE){

# Delete V char from rownames
rownames(tempDF)<-gsubfn("V","",rownames(tempDF))

# Sort data.frame according to rownames
tempDF<-tempDF[order(as.numeric(rownames(tempDF))),,drop=FALSE]

} else if(with.labels==TRUE){

rownames(labelsDF) <- labelsDF[,2]

tempDF <- as.data.frame(tempDF[rownames(labelsDF),])

}

# Check if there are any NA values and zero them before summing
tempDF[is.na(tempDF)]<-0

# Scale results from 0 to 100
tempDF.F<-(tempDF/sum(tempDF[,1]))*100*impCalcScaleF[,i]

cols=i
rows=1
for(rows in 1:nrow(tempDF.F)){
matrycaVarImp.F[rows,cols]<-tempDF.F[rows,]
}

}

# Sum row-by-row
for (rows in 1:nrow(matrycaVarImp.F)){
matrycaVarImp.F[rows,(ncol(matrycaVarImp.F)-2)]<-sum(matrycaVarImp.F[rows,])
}


# Calculate percentages
maks.sum.f <- max(matrycaVarImp.F[,(ncol(matrycaVarImp.F)-2)])

for (rows in 1:nrow(matrycaVarImp.F)){
matrycaVarImp.F[rows,(ncol(matrycaVarImp.F)-1)]<-(matrycaVarImp.F[rows,(ncol(matrycaVarImp.F)-2)])/maks.sum.f*100
}


matrycaVarImp.F <- as.data.frame(matrycaVarImp.F)

filesVarImpCols <- gsubfn(".*_VarImp_","",filesVarImp)
filesVarImpCols <- gsubfn(".txt","",filesVarImpCols)

colnames(matrycaVarImp.F)<- filesVarImpCols

names(matrycaVarImp.F)[length(matrycaVarImp.F)-2]<-"SUM"

names(matrycaVarImp.F)[length(matrycaVarImp.F)-1]<-"SUM%"

names(matrycaVarImp.F)[length(matrycaVarImp.F)]<-"ImpGrad"

# Sort data.frame according to SUM variable importance
matrycaVarImp.F<-matrycaVarImp.F[order(matrycaVarImp.F[,length(matrycaVarImp.F)-2],decreasing=TRUE),,drop=FALSE]

# # Importance gradient
for (rows in 2:nrow(matrycaVarImp.F)){

imp1 <- matrycaVarImp.F[(rows-1),(ncol(matrycaVarImp.F)-2)]
imp2 <- matrycaVarImp.F[(rows),(ncol(matrycaVarImp.F)-2)]

matrycaVarImp.F[rows,ncol(matrycaVarImp.F)] <- (imp1 - imp2)/imp1*100

}

# Add importance variables names as last column
matrycaVarImp.F$Input_no <- rownames(matrycaVarImp.F)

impCalcRes <- list("rawMeasureError"=rawMeasure, "matrixVarImp.MeasureError"=matrycaVarImp.F,
                   "model"=modelRes)

}

return(impCalcRes)

}