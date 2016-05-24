FBNormalization <- function (rawDataFileName = NULL, fishProbesFileName = NULL, normDataFileName = NULL, debugFlag = FALSE, plotFlag = FALSE, plotAndSaveFlag = FALSE)
{
	k = list.files()
	########################################
	# check the existence of the input files
	if((length(which(k == rawDataFileName)) == 0) || is.null(rawDataFileName)){
		cat("ERROR - can not find rawDataFileName: ", rawDataFileName, "\n")
		return(NULL)	
	}
 	if((length(which(k == fishProbesFileName)) == 0) || is.null(fishProbesFileName)){
		cat("ERROR - can not find fishProbesFileName:" ,  fishProbesFileName, "\n")
		return(NULL)
	}
	##############################################################################
	# check if the output file exists or not, and ask if necessary for overwriting
	while((length(which(k == normDataFileName)) != 0) || is.null(normDataFileName) || (normDataFileName == "")) {
		if(is.null(normDataFileName) || (normDataFileName == "")){
			cat("WARNING: normDataFileName is NULL.\n")
			r = "n"
		}
		else{
			cat("WARNING: ", normDataFileName, " already exists. Overwrite Yes/No/Stop?(y/n/s):")
			r = ""
			while((r != "y") && (r != "n")&& (r != "s")) {
				r = readline()
			}
		}
		if(r == "s") return()
		if(r != "y") {
			cat("Please input a new output file name: ")
			normDataFileName = readline()
		}
		if(r == "y") break
	}	
	rm(k)
	if(length(grep(".txt", normDataFileName)) == 0)
		normDataFileName = paste(normDataFileName, ".txt", sep = "")
	######################
	# read the input files
	rawData = read.table(rawDataFileName,header=TRUE,sep="\t")
	fishProbe = read.table(fishProbesFileName,header=TRUE,sep="\t")
    ##########################
    # define the output matrix
	normData = rawData
    normData[,5:dim(normData)[2]] = NA; 
	######################
    # copy raw data fields
    chrData = rawData$chrom
    nameData = names(rawData)[5:dim(rawData)[2]]
    posData = rawData$physical.position
    cytoData = as.character(rawData$cytoband)
    rawData = as.matrix(rawData[,-(1:4)])
    rowData = dim(rawData)[1]
    colData = dim(rawData)[2]
    #######################
    # copy FISH data fields
    nameFish = names(fishProbe)[6:dim(fishProbe)[2]]
    chrFish = fishProbe$chromosome
    cytoFish = fishProbe$cytoband
    cloneFish = fishProbe$BACclone
    startFish = fishProbe$start.loc
    endFish = fishProbe$end.loc
    fishProbe = as.matrix(fishProbe[,-(1:5)])
    rowFish = dim(fishProbe)[1]
    colFish = dim(fishProbe)[2]
    #########################
    # check the files headers
    flagError = FALSE
    if(is.null(chrData)){
    	flagError = TRUE
    	cat("Please name \"chrom\" the field containing the chromosomes index in the rawData file!\n")
    }
    if(is.null(posData)){
    	flagError = TRUE
    	cat("Please name \"physical.position\" the field containing the physical position of the SNP probes in the rawData file!\n")
    }
    if(is.null(cytoData)){
    	flagError = TRUE
    	cat("Please name \"cytoband\" the field containing the cytoband location of the SNP probes in the rawData file!\n")
    }
    if(is.null(chrFish)){
    	flagError = TRUE
    	cat("Please name \"chromosome\" the field containing the chromosomes index in the FISH probes file!\n")
    }
    if(is.null(cytoFish)){
    	flagError = TRUE
    	cat("Please name \"cytoband\" the field containing the cytoband location of the FISH probes in the FISH probes file!\n")
    }
    if(is.null(cloneFish)){
    	flagError = TRUE
    	cat("Please name \"BACclone\" the field containing the clone name of the FISH probes in the FISH probes file!\n")
    }
    if(is.null(startFish)){
    	flagError = TRUE
    	cat("Please name \"start.loc\" the field containing the start location of the FISH probes in the FISH probes file!\n")
    }
    if(is.null(endFish)){
    	flagError = TRUE
    	cat("Please name \"end.loc\" the field containing the start location of the FISH probes in the FISH probes file!\n")
    }
	if(flagError) return()
    ###############################################
    # test if all data have a correspondent in FISH
    # and eliminate the possible additional entries
    flagHasFish = vector(mode = "numeric", length = colData)
    newFishProbes = matrix(data = NA, nrow = rowFish, ncol = colData)
    newNameFish = vector(mode = mode(nameFish), length = colData)
    for(i in 1:colData){
        j = which(nameFish == nameData[i])
        if(length(j) >= 1){
            flagHasFish[i] = 1
            newFishProbes[,i] = fishProbe[, j]
            newNameFish[i] = nameFish[j]
        } else
            newNameFish[i] = NA
    }
    if(length(flagHasFish[flagHasFish == 0]) > 0){
        cat("WARNING: There is data with no FISH associated info: Median Centering Normalization will be used\n")
    }
    fishProbe = newFishProbes
    nameFish = newNameFish
    rm(newFishProbes, newNameFish) 
    colFish = dim(fishProbe)[2] 
    #########################################
    # Associate the CN read by the FISH probe
    # to the corresponding raw SNP value     
    snpProbe = matrix(data = NA, nrow = rowFish, ncol = colData)
    for( i in 1:rowFish){
        posProbe = posData
        posProbe[chrData != chrFish[i]] = 0
        idxStartPosProbe = which(posProbe > startFish[i])
        idxEndPosProbe = which( (posProbe < endFish[i])&(posProbe != 0) )
        if( (length(idxStartPosProbe) == 0)|(length(idxEndPosProbe) == 0) ){
            cat("WARNING: Fish probe ", cloneFish[i], " has no physical correspondence with the data...\n" )
            rm(posProbe, idxStartPosProbe, idxEndPosProbe)
            next
        }
        if(idxStartPosProbe[1] > idxEndPosProbe[length(idxEndPosProbe)]){
            idxPosProbe = idxEndPosProbe[length(idxEndPosProbe)]:idxStartPosProbe[1]
        } else {
            idxPosProbe = idxStartPosProbe[1]:idxEndPosProbe[length(idxEndPosProbe)]
        }
        allSnpProbe = rawData[idxPosProbe, ]
        if(length(idxPosProbe) == 1){
        	snpProbe[i,] = allSnpProbe
        } else{
	        for( j in 1:colData)
    	        snpProbe[i,j] = median(allSnpProbe[,j])
    	}
        rm(posProbe, idxStartPosProbe, idxEndPosProbe, idxPosProbe)
    }
    ###################################################################
    # apply k-means on each probe to determine the normalization values
    histRawData = hist(as.vector(as.matrix(rawData)), plot = FALSE, breaks = "FD")
	breaksRawData = histRawData$breaks
    normalizingSNP = matrix(data = NA, nrow = rowFish, ncol = colData)
    for(i in 1:colData){
		if(flagHasFish[i] == 0) next
		lineData = rawData[chrData != chrData[length(chrData)], i]
		clusterLineData = FBN.kmeans(inputData = lineData, minSpan = 0.2, breaksData = breaksRawData)
		for(j in 1:rowFish){
			# find the cluster that contains the snpProbe value #
			ind = clusterLineData$cluster[lineData == snpProbe[j,i]]
			normalizingSNP[j, i] = median(lineData[clusterLineData$cluster == ind[1] ] )
		}
		rm(lineData, clusterLineData)
    }
    rm(histRawData)
    ###################################
    # Check the FISH probes consistency
    for(i in 1:colData){
    	if(flagHasFish[i] == 0) next
    	cn = fishProbe[,i]
    	snp = normalizingSNP[!is.na(cn),i]
    	probe = snpProbe[!is.na(cn),i]
    	cn = cn[!is.na(cn)]
    	if(length(cn) == 1) next    	
		sortSnp = sort(snp, index.return = TRUE)
		sortCN = cn[sortSnp$ix]
		sortProbe = probe[sortSnp$ix]
		# verify if CN is in ascending order
		diffSnp = sortSnp$x[2:length(sortSnp$x)] - sortSnp$x[1:(length(sortSnp$x)-1)]
		diffCN = sortCN[2:length(sortCN)] - sortCN[1:(length(sortCN)-1)]
		diffSnp[diffSnp>0] = 1
		diffSnp[diffSnp<0] = -1
		diffCN[diffCN>0] = 1
		diffCN[diffCN<0] = -1
		flagInconsistency = FALSE
		for(j in 1:length(diffCN)){
			if(diffSnp[j] != diffCN[j])
				flagInconsistency = TRUE
		}
		if(flagInconsistency){
			flagHasFish[i] = 2
			cat("WARNING: Found inconsistent FISH info for probe: ", nameData[i], "\n")
			cat("NormalizingSNP -> fishCN -> rawSNP\n")
	    	for(j in 1:length(sortCN))
				cat(sortSnp$x[j], " -> ", sortCN[j], " -> ", sortProbe[j], "\n")
		}
	}
    
    #################################
    # Start the normalization process
    #################################
    #make three stages normalization
    # first all probes with good fish data -> flagHasFish == 1
    # second all probes with inconsistent data -> flagHasFish == 2
    # third all probes with no fish data - median centering -> flagHasFish == 0

	########################################################
    # set up the order of the CN iterations
    orderCN = c(2, 1, 3:6)
    valueCN = vector(mode = "numeric", length = length(orderCN)) #initialize the nominal values vector of the CNs
    valueCN[2] = 2
    thresholdsCN = vector(mode = "numeric", length = 6) # thresholds = {CN0-1, CN1-2, CN2-3, CN3-4, CN4-5, CN5-more}
    flagIsNormalized = vector(mode = "numeric", length = colData)
    deltaBreaksRawData = median(breaksRawData[2:length(breaksRawData)] - breaksRawData[1:(length(breaksRawData)-1)])
    breaksNormData  = vector(mode = "numeric", length = 2*max(breaksRawData)/deltaBreaksRawData)
    for(i in 1:length(breaksNormData))
    	breaksNormData[i] = (i-1)*deltaBreaksRawData
    ##############################################################################################################
    # First stage!
    # All probes with good fish data -> flagHasFish == 1
    # Normalize on the closest fishCN == 2
    for(j in 1:length(orderCN) ){
        CN = orderCN[j]
        flagNewData = FALSE
		for(i in 1:colData){
			if((flagIsNormalized[i] == 1)||(flagHasFish[i] != 1)) next
			index = which((fishProbe[,i] - CN) == 0) #look for a normalizing fish-snp pair
			if(length(index) == 0) next # current probe does not have a normalizing fish-snp pair
	    	cat("Normalize probe: ", nameData[i], "\n")
			flagNewData = TRUE
			index = index[1]
			normData[,4+i] = FBN.valueCenter(inputData = rawData[,i], normalizingValue = normalizingSNP[index,i], nominalValueCN = valueCN[CN])
			flagIsNormalized[i] = 1
			if(debugFlag){
				par(mfrow = c(2, 1))
				hist(rawData[,i], plot = TRUE, breaks = breaksNormData, main = nameData[i])
				hist(normData[,4+i], plot = TRUE, breaks = breaksNormData)
				if(plotAndSaveFlag){
					dev.print(device = bmp, width=1024, height=768, paste(nameData[i],".bmp",sep="")) 
				}
				cat("push enter to continue....")
				readline()
			}
			if(plotFlag){
				par(mfrow = c(2, 1))
				hist(rawData[,i], plot = TRUE, breaks = breaksNormData, main = nameData[i])
				hist(normData[,4+i], plot = TRUE, breaks = breaksNormData)
			}
			if(plotAndSaveFlag){
				par(mfrow = c(2, 1))
				hist(rawData[,i], plot = TRUE, breaks = breaksNormData, main = nameData[i])
				hist(normData[,4+i], plot = TRUE, breaks = breaksNormData)
				dev.print(device = bmp, width=1024, height=768, paste(nameData[i],".bmp",sep="")) 
			}
		}
		if(!flagNewData) next
		cat("\nK-MEANS for group CN", CN, "\n")
        allNormData = as.matrix(normData[chrData != chrData[length(chrData)], 4+which(flagIsNormalized == 1)])
        allNormData = as.numeric(as.vector(allNormData))
        if(debugFlag | plotFlag | plotAndSaveFlag){
        	par(mfrow = c(1, 1))
        	histAllNormData = hist(allNormData, plot = TRUE, breaks = breaksNormData)
        }
        #set the minSpan to 0.3 for the final clustering...
        clusterAllNormData = FBN.kmeans(inputData = allNormData, minSpan = 0.3, breaksData = NULL)
		#plot(allNormData, col = clusterAllNormData$cluster)
		tempNominalCN = abs(clusterAllNormData$centers - 2)
		idxCN2 = which(tempNominalCN == min(tempNominalCN))
		# merge all clusters that are inferior to CN2 (due to possible partial hybridizations)
		clusterAllNormData$cluster[clusterAllNormData$cluster < idxCN2] = idxCN2 - 1
		indexCNs = c()
		for(i in 1:6){
			if((idxCN2-2+i) < 1) next
			if((idxCN2-2+i) > length(clusterAllNormData$centers))
				break
			indexCNs = c(indexCNs, idxCN2-2+i)
			#valueCN[i] = clusterAllNormData$centers[idxCN2-2+i]
			valueCN[i] = median(allNormData[clusterAllNormData$cluster == indexCNs[i]])
			if(i == 1)
				thresholdsCN[i] = min(allNormData[clusterAllNormData$cluster == (idxCN2-2+i)])
			else
				thresholdsCN[i] = (min(allNormData[clusterAllNormData$cluster == (idxCN2-2+i)]) + max(allNormData[clusterAllNormData$cluster == (idxCN2-2+i-1)]))/2
		}
		if(length(which(flagIsNormalized[flagHasFish == 1] == 0)) > 0)
			rm(allNormData, clusterAllNormData)
		if(debugFlag | plotFlag | plotAndSaveFlag)
			rm(histAllNormData)
		cat("\nCN values: ", valueCN, "\n\n")
        if(debugFlag){
        	cat("push enter to continue....")
			readline()
		}
    }
    ##############################################################################################################
    # Second stage!
    # Normalize on each available fish data and verify all others. 
    # Keep the normalization with the highest success score
    for(i in 1:colData){
    	if((flagIsNormalized[i] == 1)||(flagHasFish[i] != 2)) next
    	score = vector(mode = "numeric", length = rowFish)
    	cat("Normalize probe: ", nameData[i], "\n")
    	#normalize on each fish info
    	for(j in 1:rowFish){
    		if(is.na(fishProbe[j, i])) next
			# need to normalize only the normalizingSNP to calculate the score. Then normalize the data for a maximum score...
			tmpNormSNP = FBN.valueCenter(inputData = normalizingSNP[,i], normalizingValue = normalizingSNP[j,i], nominalValueCN = valueCN[fishProbe[j,i]])
			#calculate score for current normalization
			tmpCN = vector(mode = "numeric", length = rowFish)
			score[j] = 0
			for(k in 1:rowFish){
				if(is.na(fishProbe[k, i])) next
				tmp = thresholdsCN - tmpNormSNP[k]
				idx1 = which(tmp <= 0)
				idx2 = which(tmp >= 0)
				if(length(idx2)>0)
					tmpCN[k] = idx2[1]-1
				else 
					tmpCN[k] = 5
				if(tmpCN[k] == fishProbe[k, i])
					score[j] = score[j] + 1
			}
		}
		idxNorm = which(score == max(score))
		cat(" with success scores", score, "\n")
		normData[,4+i] = FBN.valueCenter(inputData = rawData[,i], normalizingValue = normalizingSNP[idxNorm[1],i], nominalValueCN = valueCN[fishProbe[idxNorm[1],i]])
		flagIsNormalized[i] = 1
		if(debugFlag){
			par(mfrow = c(2, 1))
			hist(rawData[,i], plot = TRUE, breaks = breaksNormData, main = nameData[i])
			hist(normData[,4+i], plot = TRUE, breaks = breaksNormData)
			#dev.print(device = bmp, width=1024, height=768, paste(nameData[i],".bmp",sep="")) 
			cat("push enter to continue....")
			readline()
		}		
		if(plotFlag){
			par(mfrow = c(2, 1))
			hist(rawData[,i], plot = TRUE, breaks = breaksNormData, main = nameData[i])
			hist(normData[,4+i], plot = TRUE, breaks = breaksNormData)
		}
		if(plotAndSaveFlag){
			par(mfrow = c(2, 1))
			hist(rawData[,i], plot = TRUE, breaks = breaksNormData, main = nameData[i])
			hist(normData[,4+i], plot = TRUE, breaks = breaksNormData)
			dev.print(device = bmp, width=1024, height=768, paste(nameData[i],".bmp",sep="")) 
		}
    }
    ##############################################################################################################
    # Third stage!
    # Classical Median Centering Normalization - The median value of the raw data is normalized in CN2
    for(i in 1:colData){
    	if((flagIsNormalized[i] == 1)||(flagHasFish[i] != 0)) next
    	cat("Normalize probe: ", nameData[i], "\n")
    	normData[,4+i] = FBN.valueCenter(inputData = rawData[,i], normalizingValue = median(rawData[chrData != chrData[length(chrData)],i]), nominalValueCN = 2)
		flagIsNormalized[i] = 1
    }
    ############################################################
    # Make the final clusterization and determine the thresholds
	cat("\nFinal K-MEANS for all data", "\n")
    allNormData = as.matrix(normData[chrData != chrData[length(chrData)], 4+which(flagIsNormalized == 1)])
    allNormData = as.numeric(as.vector(allNormData))
	if(plotFlag | plotAndSaveFlag){
    	par(mfrow = c(1, 1))
    	histAllNormData = hist(allNormData, plot = TRUE, breaks = breaksNormData)
    }
	clusterAllNormData = FBN.kmeans(inputData = allNormData, minSpan = 0.3, breaksData = NULL)
	#plot(allNormData, col = clusterAllNormData$cluster)
	tempNominalCN = abs(clusterAllNormData$centers - 2)
	idxCN2 = which(tempNominalCN == min(tempNominalCN))
	######################################################################################
	# merge all clusters that are inferior to CN2 (due to possible partial hybridizations)
	clusterAllNormData$cluster[clusterAllNormData$cluster < idxCN2] = idxCN2 - 1
	indexCNs = c()
	for(i in 1:6){
		if((idxCN2-2+i) < 1) next
		if((idxCN2-2+i) > length(clusterAllNormData$centers))
			break
		indexCNs = c(indexCNs, idxCN2-2+i)
		valueCN[i] = median(allNormData[clusterAllNormData$cluster == indexCNs[i]])
		if(i == 1)
			thresholdsCN[i] = min(allNormData[clusterAllNormData$cluster == (idxCN2-2+i)])
		else
			thresholdsCN[i] = (min(allNormData[clusterAllNormData$cluster == (idxCN2-2+i)]) + max(allNormData[clusterAllNormData$cluster == (idxCN2-2+i-1)]))/2
	}
	cat("\nCN values: ", valueCN, "\n\n")
    ############################################################
    # determine the thresholds
    # fit gaussian on each cluster and find their meeting points
    if(length(indexCNs) == 6)
    	indexCNs = indexCNs[1:5]
    meanCNs = vector(mode = "numeric", length = length(indexCNs))
    stdCNs = vector(mode = "numeric", length = length(indexCNs))
    for(i in 1:length(indexCNs)){
	    meanCNs[i] = mean(allNormData[clusterAllNormData$cluster == indexCNs[i]])
		stdCNs[i] = sd(allNormData[clusterAllNormData$cluster == indexCNs[i]])
    }
    #################################################################
    # thresholds = {CN0-1, CN1-2, CN2-3, CN3-4, CN4-5, CN5-6 or more}
    thresholds = vector(mode = "numeric", length = (length(indexCNs)+1))
    thresholds[1] = meanCNs[1] - 2*stdCNs[1]
    thresholds[length(indexCNs)+1] = meanCNs[length(indexCNs)] + 2*stdCNs[length(indexCNs)]
    for(i in 2:(length(indexCNs))){
    	thresholds[i] = meanCNs[i] - (meanCNs[i] - meanCNs[i-1])*stdCNs[i]/(stdCNs[i] + stdCNs[i-1])
    }
    for(i in 1:length(thresholds)){
	    cat("threshold CN", i-1, " to CN", i, " = ", thresholds[i], "\n")
	}
	write.table(normData, file = normDataFileName, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
	thresholdsFileName = strsplit(normDataFileName, ".txt");
	thresholdsFileName = paste(thresholdsFileName, "_thresholds.txt", sep = "")
	b = c("CN1", "CN2", "CN3", "CN4", "CN5", "")
	b = c(b[1:length(indexCNs)], "")
	a = valueCN
	a = c(a[1:length(indexCNs)], "")

	bb = c("ThrCN0-1", "ThrCN1-2", "ThrCN2-3", "ThrCN3-4", "ThrCN4-5", "ThrCN5-6")
	bb = bb[1:(length(indexCNs)+1)]
	bb[length(bb)] = paste(bb[length(bb)], " or more", sep = "")
	aa = thresholds
	aa = aa[1:(length(indexCNs)+1)]
	outInfo = rbind(b, a, bb, aa)
	
	write.table(outInfo, file = thresholdsFileName, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
	cat("FBNormalization done!\n")
	return()
}
