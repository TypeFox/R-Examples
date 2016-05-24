# TODO: Add comment
# 
# Author: ahrnee-adm
###############################################################################


### INIT
if(!grepl("SafeQuant\\.Rcheck",getwd())){
	setwd(dirname(sys.frame(1)$ofile))
}
source("initTestSession.R")
### INIT END

### TEST FUNCTIONS

#install.packages("/Users/erikahrne/dev/R/workspace/SafeQuant/", repos = NULL, type="source")
#library(SafeQuant)


### INIT
skylineExportFile <- "testData/skyline_dilution_curve.csv"

### parse
skylineData <- read.csv(skylineExportFile,sep=",")

### HACK replace concentration 
#if(F){
	#concentration <- gsub("B15\\-01....\\_NRX..\\_WB_GFPl:*_*1:" ,"",as.character(skylineData$Replicate.Name))
	concentration <- gsub("B15\\-0...._.*:" ,"",as.character(skylineData$Replicate.Name))
	#concentration <- gsub("B15\\-[0-9]*" ,"",as.character(skylineData$Replicate.Name))
	#concentration <- gsub("^_" ,"",concentration)
	concentration <- gsub("_WB.*","",concentration)
	concentration <- gsub("B15.*","",concentration)
	concentration <- gsub("_.*","",concentration)
	concentration <- 1/as.numeric(concentration)
	concentration[is.na(concentration)] <- 0
		
	#repConc <- data.frame( rev(c(79.55384267,16.04907907,3.352755873,0.71222776,0.278683489,0)),row.names=sort(unique(concentration)))
	repConc <- data.frame( rev(c(68.80773399,	14.13167898,	 2.878151306,	 0.608765369,	 0.167507442, 0)),row.names=sort(unique(concentration)))

	concentration <- repConc[as.character(concentration),]
	rev(sort(unique(concentration)))
	data.frame(concentration,skylineData$Replicate.Name)
	
#}
	
skylineData <- cbind(skylineData,concentration)
skylineData <- skylineData[!grepl("N",skylineData$Area),]
skylineData$Area <- as.numeric(as.character(skylineData$Area))

### rollUp -> create ExpressionSet
nbReplicates <- 3
peptideChargeMzStateConc <-  paste(skylineData$Peptide.Sequence,skylineData$Precursor.Charge,round(skylineData$Precursor.Mz),skylineData$concentration,sep="_")
uniquePeptideChargeStateConc <-  unique(peptideChargeMzStateConc)

expressionMatrix <- data.frame(matrix(nrow=length(unique(unique(peptideChargeMzStateConc))),ncol=nbReplicates),row.names=unique(peptideChargeMzStateConc))
featureAnnotations <- data.frame()
for(pcsc in uniquePeptideChargeStateConc){
	
	peptideChargeMzStateConcSubset <- skylineData[peptideChargeMzStateConc %in% pcsc,]
	
	### add feature data
	c <- 1
	for(rep in unique(peptideChargeMzStateConcSubset$Replicate.Name)){
		expressionMatrix[pcsc,c] <-sum(peptideChargeMzStateConcSubset[peptideChargeMzStateConcSubset$Replicate.Name %in% rep,]$Area)
			c <- c+1
	}

	### add assay data
	featureAnnotations  <- rbind(featureAnnotations,peptideChargeMzStateConcSubset[1,c(1,2,4,5,9,13)])
	
}

featureAnnotations <-cbind(featureAnnotations,dilutionCurveId=paste(featureAnnotations$Peptide.Sequence,featureAnnotations$Precursor.Charge,round(featureAnnotations$Precursor.Mz),sep="_")) 
expressionMatrix[expressionMatrix == 0] <- NA 
row.names(featureAnnotations) <- unique(peptideChargeMzStateConc)	
expDesign <- data.frame(condition=rep(1,3),isControl=rep(T,3))
names(expressionMatrix) <- rownames(expDesign)

### rollUp -> create ExpressionSet END
esetCalibCurve <- createExpressionDataset(expressionMatrix=as.matrix(expressionMatrix),expDesign=expDesign,featureAnnotations=featureAnnotations)
### INIT END


testCalibrationCurve <- function(){
	
	cat(" --- testcalibrationCurve --- \n")
	
	calibCurve <- calibrationCurve(esetCalibCurve[fData(esetCalibCurve)$dilutionCurveId == unique(fData(esetCalibCurve)$dilutionCurveId)[2], ])
	stopifnot(all.equal(round(calibCurve$lod,2),0.44))
	cat(" --- testcalibrationCurve: PASS ALL TEST  --- \n")
	
}

### RUN TESTS

testCalibrationCurve()

### RUN TESTS END

### GRAPHICS

#indices <- as.character(fData(esetCalibCurve)$dilutionCurveId)
#plot(calibrationCurve(esetCalibCurve[indices == unique(indices)[2], ],method="low"), xlab="Concentration (fmol/ul)")
#plot(calibrationCurve(esetCalibCurve[indices == unique(indices)[2], ],method="blank"), xlab="Concentration (fmol/ul)")

print("DONE")

