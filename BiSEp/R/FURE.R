# Functional Redundancy Algorithm Build V2.0
#
# Date of V1.0: November - December 2012
# Date of V2.0: June 2013	
# Author: M Wappett
# Decription: Input a matrix of continuous data with genes as rows and samples as columns.  Algorithm may take some considerable time to run (~ 10 hours for 80,000 gene pairs as results)

# Read command line arguments
FURE <- function(data=data, inputType=inputType)
{

# Define confidence criteria
if(missing(data)) stop("Need to input expression data matrix")
if(missing(inputType)) stop("Need to specify sample type")

# Is the input matrix from the SlinG or BEEM algorithms?
colNum <- 0
if(inputType == "BIGEE")
{
	colNum <- 4
}
else if(inputType == "BEEM")
{
	colNum <- 11
}

# Do GO term functional redundancy mapping
tab1 <- merge(toTable(org.Hs.eg.db::org.Hs.egSYMBOL2EG), toTable(org.Hs.eg.db::org.Hs.egGO))
w1 <- which(data[,1] %in% tab1[,2])
x8 <- data[w1,]
w2 <- which(x8[,2] %in% tab1[,2])
x9 <- x8[w2,]
unG4 <- unique(c(x9[,1], x9[,2]))
tab1 <- merge(toTable(org.Hs.eg.db::org.Hs.egSYMBOL2EG[unG4]), toTable(org.Hs.eg.db::org.Hs.egGO))
tab2 <- select(GO.db::GO.db, tab1$go_id, "TERM", "GOID")
tab3 <- cbind(tab1, tab2)
tab4 <- unique(tab1[,1:2])
redundantIDs <- character(0)
redundantTerms <- character(0)
outList <- list(0)
for(i in 1:dim(x9)[1])
{
	s1 <- subset(tab3, tab3[,2] == x9[i,1])
	s2 <- subset(tab3, tab3[,2] == x9[i,2])
	outList[[i]] <- c(s1[1,1], s2[1,1])
	m1 <- merge(s1, s2, by.x="TERM", by.y="TERM")
	if(dim(m1)[1] > 0)
	{
		redundantIDs[i] <- paste(m1[,4], collapse = " / ")
		redundantTerms[i] <- paste(m1[,1], collapse = " / ")
	}
	else
	{
		redundantIDs[i] <- "No Redundancy"
		redundantTerms[i] <- "No Redundancy"
	}
}
x9$redundant_ids <- redundantIDs
x9$redundant_terms <- redundantTerms

unG5 <- as.list(tab4[,1])
# Perform semantic similarity mapping from gene Ontologies (this may take some considerable time - up to 8 hours; suggest running overnight)
v1 <- mgeneSim(unG5, ont="MF", organism="human")
v2 <- mgeneSim(unG5, ont="BP", organism="human")
v3 <- mgeneSim(unG5, ont="CC", organism="human")

mfScore <- character(0)
for(i in 1:dim(x9)[1])
{
	w1 <- which(tab4[,2] %in% x9[i,1])
	w2 <- which(tab4[,2] %in% x9[i,2])
	w3 <- tab4[w1,1]
	w4 <- tab4[w2,1]
	w5 <- which(colnames(v1) == w4 | colnames(v1) == w3 )
	w6 <- which(rownames(v1) == w4 | rownames(v1) == w3 )
	w7 <- length(c(w5, w6))
	if(w7 == 4)
	{
		mfScore[i] <- v1[w5[1],w6[2]]
	}
	else
	{
		mfScore[i] <- "No Score"
	}
}

bpScore <- character(0)
for(i in 1:dim(x9)[1])
{
	w1 <- which(tab4[,2] %in% x9[i,1])
	w2 <- which(tab4[,2] %in% x9[i,2])
	w3 <- tab4[w1,1]
	w4 <- tab4[w2,1]
	w5 <- which(colnames(v2) == w4 | colnames(v2) == w3 )
	w6 <- which(rownames(v2) == w4 | rownames(v2) == w3 )
	w7 <- length(c(w5, w6))
	if(w7 == 4)
	{
		bpScore[i] <- v2[w5[1],w6[2]]
	}
	else
	{
		bpScore[i] <- "No Score"
	}
}

ccScore <- character(0)
for(i in 1:dim(x9)[1])
{
	w1 <- which(tab4[,2] %in% x9[i,1])
	w2 <- which(tab4[,2] %in% x9[i,2])
	w3 <- tab4[w1,1]
	w4 <- tab4[w2,1]
	w5 <- which(colnames(v3) == w4 | colnames(v3) == w3 )
	w6 <- which(rownames(v3) == w4 | rownames(v3) == w3 )
	w7 <- length(c(w5, w6))
	if(w7 == 4)
	{
		ccScore[i] <- v3[w5[1],w6[2]]
	}
	else
	{
		ccScore[i] <- "No Score"
	}
}

x9$MolecularFunctionScore <- mfScore
x9$BiologicalProcessScore <- bpScore
x9$CellularComponentScore <- ccScore
w1 <- which(x9[,colNum] == "No Redundancy")
x10 <- x9[-w1,]
outList <- list(allPairsScored=x9, funcRedundantPairs=x10)
return(outList)
}
