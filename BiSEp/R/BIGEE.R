# Bimodal Gene Expression Exclusivity Build V3.0
#
# Date of V1.0: November - December 2012
# Date of V2.0: April 2013
# Date of V3.0: May - June 2013	
# Author: M Wappett
# Decription: Takes as input the list object output from the BISEP function

# Load relevant libraries
BIGEE <- function(bisepData=data, sampleType=c("cell_line", "cell_line_low", "patient", "patient_low"))
{

# Define confidence criteria
if(missing(bisepData)) stop("Need to input expression data matrix")
if(missing(sampleType)) stop("Need to specify sample type")
if(sampleType == "cell_line")
{
	print("Selected CELL LINE sample type")
	pI <- 0.5
	dTA <- 2.5
	bI <- 0.7
	numALL <- 1.9
}
else if(sampleType == "patient")
{	
	print("Selected PATIENT sample type")
	pI <- 0.5
	dTA <- 2.5
	bI <- 0.5
	numALL <- 1.9
}
else if(sampleType == "cell_line_low")
{	
	print("Selected CELL LINE sample type")
	pI <- 0.5
	dTA <- 3.5
	bI <- 1.1
	numALL <- 1.9
}
else if(sampleType == "patient_low")
{	
	print("Selected PATIENT sample type")
	pI <- 0.5
	dTA <- 3
	bI <- 0.9
	numALL <- 1.9
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

print("Subsetting bimodal index")

subBiIndex <- subset(biIndex, biIndex[,6] > bI & biIndex[,5] < pI & biIndex[,4] > dTA)
subBiIndex2 <- subBiIndex[order(-subBiIndex[,6]),]

print("Filtering")

# Calculate the mid-points for all genes identified as being bimodal
midpoint <- big.model[,1]
names(midpoint) <- rownames(data2)
w1 <- rownames(subBiIndex)
if(length(w1) <= 1)
{
	stop("No bimodal genes in input matrix")
}
midpointBI <- midpoint[w1]

data3 <- data2[w1,]
rownames(data3) <- w1

med1 <- apply(data3, 1, function(x) mean(as.numeric(x[1:dim(data3)[2]])))
w1 <- which(med1 > 2)
data3 <- data3[w1,]
w1 <- rownames(data3)
midpointBI <- midpointBI[w1]

print("Setting up synthetic lethal detection")

# Set data up for pure synthetic lethal detection
list3 <- as.list(as.data.frame(t(data3)))
outList <- list(data.frame=(0))
ind1 <- 0
names1 <- names(list3)

print("Running SL detection")
# Run pure synthetic lethal detection
subBiIndex3 <- subBiIndex[rownames(data3),]
rangeTab <- rbind(c(0.55, 0.5), c(0.5, 0), c(0.45, 0.5), c(0.4, 1), c(0.35, 1.5), c(0.30, 2), c(0.25, 2.5), c(0.20, 3), c(0.15, 3.5), c(0.10, 4), c(0.05, 4.5), c(0.00, 5))
print(paste("Number of bimodal genes: ", length(list3), sep=""))
pb <- txtProgressBar(min=0, max=length(list3), style=3)
for(i in 1:length(list3))
{
	range1 <- range(list3[[i]])
	range1_2 <- ((range1[2] - range1[1])/100)* (rangeTab[which.min(abs(rangeTab[,1] -subBiIndex3[i,5])),2])
	midpoint1 <- midpointBI[i]
	for(j in 1:(length(list3)))
	{
		range2 <- range(list3[[j]])
		range2_2 <- ((range2[2] - range2[1])/100)* (rangeTab[which.min(abs(rangeTab[,1] -subBiIndex3[j,5])),2])
		midpoint2 <- midpointBI[j]
		midpoint3 <- midpoint1 - range1_2
		midpoint4 <- midpoint2 - range2_2
		num1 <- which(list3[[i]] < midpoint3)
		pair2 <- list3[[j]][num1]
		num2 <- which(pair2 < midpoint4)
		num3 <- which(list3[[j]] < midpoint4)
		num4 <- (sum(length(num1), length(num3))/(length(list3[[1]])*2)) * 100
		len1 <- c(length(num1), length(num3))
		perLL <- length(num2)/length(list3[[i]])*100
		sd1 <- sd(len1)
		rank1 <- (num4 / sd1)
		if(perLL == 0 & num4 > numALL & num4 < 50)
		{
			ind1 <- ind1 + 1
			BBMM <- c(names1[i], names1[j])
			outList[[ind1]] <- BBMM
		}
		else
		{
			#Do Nothing
		}		
	}
	Sys.sleep(0.5)
	setTxtProgressBar(pb, i)
	
}
Sys.sleep(1)
close(pb)
print("Summarising...")
x2 <- do.call("rbind", outList)
x3 <- as.data.frame(x2[,1:2], stringsAsFactors=FALSE)
x4 <- x3[x3[,1] <= x3[,2],]
colnames(x4) <- c("gene", "gene2")

biIndex$gene <- rownames(biIndex)
colnames(biIndex)[7] <- "gene"
x5 <- merge(x4, biIndex, by.x="gene", by.y="gene")
colnames(biIndex)[7] <- "gene2"
x6 <- merge(x5, biIndex, by.x="gene2", by.y="gene2")
x7 <- x6[,1:2]
score1 <- ((x6[,8] + x6[,14]) * (x6[,7] + x6[,13]))* (x6[,6] + x6[,12])
x7$score <- score1
x7 <- x7[order(-x7[,3]),]
return(x7)
print("Complete")
}
