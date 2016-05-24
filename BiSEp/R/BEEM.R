# Mutation and Gene Expression Automated Detected Algorithm Build V2.0
#
# Date of V1.0: February - March 2013
# Date of V2.0: May - June 2013	
# Author: M Wappett
# Decription: Input a matrix of continuous data with genes as rows and samples as columns.  Input a matrix of discreet mutation data ("WT" or "MUT" calls) with genes as rows and samples as columns.  Algorithm may take some considerable time to run (1 hour for 20,500 genes across 800 samples).

BEEM <- function(bisepData=data, mutData=mutData, sampleType=c("cell_line", "cell_line_low", "patient", "patient_low"), minMut=10)
{
# Define confidence criteria
if(missing(bisepData)) stop("Need to input BISEP object")
if(missing(mutData)) stop("Need to input mutation data matrix")
if(missing(sampleType)) stop("Need to specify sample type")
print(paste("Minimum number of mutations considered for each gene is: ", minMut, sep=""))
if(sampleType == "cell_line")
{
	print("Selected CELL LINE sample type")
	pI <- 0.5
	dTA <- 2.5
	bI <- 0.7
}
else if(sampleType == "patient")
{	
	print("Selected PATIENT sample type")
	pI <- 0.5
	dTA <- 2.5
	bI <- 0.5
}
else if(sampleType == "cell_line_low")
{	
	print("Selected CELL LINE LOW sample type")
	pI <- 0.5
	dTA <- 3.5
	bI <- 1.1
}
else if(sampleType == "patient_low")
{	
	print("Selected PATIENT LOW sample type")
	pI <- 0.5
	dTA <- 3
	bI <- 0.9
}
else
{
	stop("Don't recognise sample type - please review options")
}

# Extract objects from input + stop if object incorrect

if("BISEP" %in% names(bisepData) && "BI" %in% names(bisepData) && "DATA" %in% names(bisepData))
{
	biIndex <- bisepData$BI
	big.model <- bisepData$BISEP
	data2 <- bisepData$DATA
}
else
{
	stop("Input object isn't from BISEP function")
}

subBiIndex <- subset(biIndex, biIndex[,6] > bI & biIndex[,5] < pI & biIndex[,4] > dTA)
subBiIndex2 <- subBiIndex[order(-subBiIndex[,6]),]

midpoint <- big.model[,1]
names(midpoint) <- rownames(data2)
w1 <- rownames(subBiIndex)
midpointBI <- midpoint[w1]

if(length(w1) <= 1)
{
	stop("No bimodal genes in input matrix")
}

data3 <- data2[w1,]
rownames(data3) <- w1

med1 <- apply(data3, 1, function(x) mean(as.numeric(x[1:dim(data3)[2]])))
w1 <- which(med1 > 2)
data3 <- data3[w1,]
w1 <- rownames(data3)
midpointBI <- midpointBI[w1]

# Perform analysis of mutation data, and then integrate with the gene expression analysis outcomes
mutData <- as.data.frame(mutData, stringsAsFactors=FALSE)

# Check size of objects + stop if both are empty
if(dim(mutData)[1] == 0)
{
	stop("No genes with mutation frequency higher than that specified")
}
if(dim(data3)[1] == 0)
{
	stop("No bimodal genes in input matrix")
}

# Evaluate the significance of this
w1 <- which(colnames(mutData) %in% colnames(data3))
mutData2 <- mutData[,w1]
c1 <- apply(mutData2, 1, function(x) length( x [x == "MUT"] ))
mutData2$mutCount <- c1
mutData3 <- subset(mutData2, mutData2$mutCount >= minMut)
w1 <- which(colnames(data3) %in% colnames(mutData3))
if(length(w1) == 0) { stop("No overlapping sample names") }
data4 <- data3[,w1]

# Check size of objects + stop if both are empty
if(dim(mutData3)[1] == 0)
{
	stop("No genes with mutation frequency higher than that specified")
}
if(dim(data3)[1] == 0)
{
	stop("No bimodal genes in input matrix")
}

# Order mutation data columns the same as expression data columns.
mutData3 <- mutData3[,colnames(data4)]
print(paste("Number of bimodal expression genes : ", dim(data4)[1], sep=""))
print(paste("Number of mutation genes wih frequency greater than", minMut, ":", dim(mutData3)[1], sep=" "))
list1 <- as.list(as.data.frame(t(data4)))
names1 <- colnames(data4)
names2 <- colnames(mutData3)
outList <- list(data.frame(0))
ind1 <- 0
pb <- txtProgressBar(min=0, max=length(list1), style=3)
for(i in 1:length(list1))
{
	data5 <- list1[[i]][order(list1[[i]])]
	names1_1 <- names1[order(list1[[i]])]
	mutData4 <- mutData3[,names1_1]
	mpNum <- midpointBI[[i]] - (5*((max(data5) - min(data5))/100))
	sampleMidIndex <- which.min(abs(data5 - mpNum))
	lower <- data5[1:sampleMidIndex]
	upper <- data5[(sampleMidIndex+1): length(data5)]
	contTab <- matrix(ncol=2, nrow=2)
	valVec <- numeric(0)
	list2 <- as.list(as.data.frame(t(mutData4)))
	for(j in 1:length(list2))
	{
		lowerM <- list2[[j]][1:sampleMidIndex]
		upperM <- list2[[j]][(sampleMidIndex+1):length(list2[[j]])]
		count1 <- length(which(lowerM == "MUT"))
		count2 <- length(which(upperM == "MUT"))
		count3 <- count1 + count2
		contTab[1,1] <- count1
		contTab[1,2] <- count2
		contTab[2,1] <- length(lower) - count1
		contTab[2,2] <- length(upper) - count2
		per1 <- (count1/length(lower))*100
		per2 <- (count2/length(upper))*100
		per3 <- per1 - per2
		pV <- fisher.test(contTab)
		ind1 <- ind1 + 1
		if(per3 > 0)
		{
			outList[[ind1]] <- c(names(list1)[i], rownames(mutData4)[j], count1, count2, pV[[1]], per1, per2, length(lowerM), length(upperM), "LowEnriched")
		}
		else
		{
			outList[[ind1]] <- c(names(list1)[i], rownames(mutData4)[j], count1, count2, pV[[1]], per1, per2, length(lowerM), length(upperM), "HighEnriched")
		}
	}
	Sys.sleep(0.5)
	setTxtProgressBar(pb, i)
}
Sys.sleep(1)
close(pb)
print("Summarising...")
genePairs2 <- do.call("rbind", outList)
genePairs2 <- as.data.frame(genePairs2, stringsAsFactors=FALSE)
genePairs2 <- genePairs2[order(genePairs2[,5]),]
genePairs3 <- subset(genePairs2, genePairs2[,10] == "HighEnriched" & genePairs2[,5] < 0.25)
genePairs3 <- as.data.frame(genePairs3, stringsAsFactors=FALSE)
genePairs3[,5] <- as.numeric(genePairs3[,5])
genePairs4 <- genePairs3[order(genePairs3[,5]),]
colnames(genePairs4) <- c("Gene1", "Gene2", "LowerExpressionMutationCount", "HighExpressionMutationCount", "Fishers P Value", "Percentage of lower samples mutated", "Percenage of high samples mutated", "Size of low expression population", "Size of high expression population", "Enrichment Status")
return(genePairs4)
print("Complete")
}
