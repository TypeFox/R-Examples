
dataPreprocess <- function(trainMatryca_nr, testMatryca_nr, labelsFrame, lk_col, lk_row, with.labels){

orig_input_no <- data.frame("Orig_labels"=c(1:(lk_col-1)))

lk_col_at_start <- lk_col

xTrain <- as.data.frame(trainMatryca_nr[,1:(ncol(trainMatryca_nr)-1)])
yTrain <- as.data.frame(trainMatryca_nr[,ncol(trainMatryca_nr)])

xTest <- as.data.frame(testMatryca_nr[,1:(ncol(testMatryca_nr)-1)])
yTest <- as.data.frame(testMatryca_nr[,ncol(trainMatryca_nr)])

trainMatryca_nr <- xTrain
testMatryca_nr <- xTest

preprocessList<-list()
# A) Check for near zero variance predictors 

# Flag as near zero if:
# 1. the percentage of unique values is less than 20% and
# 2. the ratio of the most frequent to the second most frequent value is greater than 20,

nearZeroVarMyInp <- nearZeroVar(trainMatryca_nr)

if(length(nearZeroVarMyInp)>=1){

outfile_nearZeroVarInp<-data.frame(nearZeroVarMyInp)

# Remove near zero variance columns
labelsFrame <- as.data.frame(labelsFrame[setdiff(rownames(labelsFrame),nearZeroVarMyInp),])
orig_input_no <- as.data.frame(orig_input_no[setdiff(rownames(orig_input_no),nearZeroVarMyInp),])

trainMatryca_nr <- subset(trainMatryca_nr, select=-nearZeroVarMyInp)
testMatryca_nr <- subset(testMatryca_nr, select=-nearZeroVarMyInp)


lk_col = ncol(trainMatryca_nr)

outfile_nearZeroVarData<-data.frame(trainMatryca_nr)

lk_col = ncol(trainMatryca_nr)
lk_row = nrow(trainMatryca_nr)

}

# B) Check for susceptibility to multicollinearity

# Calculate correlation matrix
trainMatryca_nr_Corr <- cor(trainMatryca_nr)

# Find variables with correlation 0.9 or more
highCorr <- findCorrelation(trainMatryca_nr_Corr, 0.90)

if(length(highCorr)>=1){

outfile_highCorrInp<-data.frame(highCorr)


# Remove columns from list highCorr
labelsFrame <- as.data.frame(labelsFrame[setdiff(rownames(labelsFrame),highCorr),])
orig_input_no <- as.data.frame(orig_input_no[setdiff(rownames(orig_input_no),highCorr),])

transposedLabelsFrame <- t(labelsFrame)
trainMatryca_nr <- subset(trainMatryca_nr, select=-highCorr)
testMatryca_nr <- subset(testMatryca_nr, select=-highCorr)

lk_col = ncol(trainMatryca_nr)


trainMatryca_nr <- data.frame(trainMatryca_nr)
colnames(trainMatryca_nr) <- labelsFrame[,1]
testMatryca_nr <- data.frame(testMatryca_nr)
colnames(testMatryca_nr) <- labelsFrame[,1]


lk_col = ncol(trainMatryca_nr)
lk_row = nrow(trainMatryca_nr)

}

exportlabelsFrame <- data.frame(orig_input_no,labelsFrame[,1])
colnames(exportlabelsFrame)<-c("Orig Input No","Labels")

trainMatryca_nr <- as.data.frame(cbind(trainMatryca_nr,yTrain))
testMatryca_nr <- as.data.frame(cbind(testMatryca_nr,yTest))

preprocessList <- list("trainMatryca"=trainMatryca_nr, "testMatryca"=testMatryca_nr, "labelsDF"=exportlabelsFrame)
  
return(preprocessList) 

}