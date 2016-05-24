###### This file contains functions used in 'meta3dMainF.R'
###### Author: Gang Wu (wggucas@gmail.com)
###### Lab: John Hogenesch's lab (http://hogeneschlab.org/)
######======================================================================================================================================
##extract time information from given columns
extractTimeColF <- function(design_timeColm, notekeys=c("design_hrColm", "hour"), designD, IDS) 
{
    noteA <- c("The setted column number of '", notekeys[1], "' is larger ",
               "than the total column number of designfile.",
               "The total column number of designfile is ")
    noteB <- c(". The reason may be mistakenly setting ",
               "this parameter, or mistakenly setting 'filestyle'.\n")

	if ( is.numeric(design_timeColm) & (length(design_timeColm) == 1) & (design_timeColm > 0) )
	{
		if ( design_timeColm > ncol(designD) )
		{
			stop(c(noteA, ncol(designD), noteB) )
		}
		return(designD[IDS, design_timeColm])
	}  else  {
		stop(c("Please set '", notekeys[1], "' as a positive numeric value ",
		   "corresponding to the column number storing timepoints-", notekeys[2],
		   " information in the 'designfile'.\n") )
	}
}
######---------------------------------------
##extract numeric number from strings
getNumber <- function(strvector)
{
    ###extract numeric number from strings
    strvector <- as.character(strvector)
    strL <- strsplit(strvector, "")
    numvector <- sapply(strL, function(z) {
                            zindex <- grep("\\d", z)
                            if (length(zindex) > 1) {
                                zout <- paste(z[zindex[1]:zindex[length(zindex)]], collapse="")
                                return(as.numeric(zout))
                            }   else if (length(zindex) == 1)  {
                                return(as.numeric(z[zindex]))
                            }   else {
                                return(NA)
                            }
                        })
    return(numvector)
}
######---------------------------------------
##check the column order as positive numeric value
checkColmF <- function(design_Colm, notekeys="")
{
	if ( is.numeric(design_Colm) & (length(design_Colm) == 1) & (design_Colm > 0) )
    {
        return(TRUE)
    }  else  {
        stop(c("Please set '", notekeys, "' as a positive numeric value ",
               "corresponding to the column number storing libraryID in ",
               "the 'designfile'.\n") )
    }
}
######---------------------------------------
##check the column order is no larger than the number of columns of 'designD'
checkColmF2 <- function(design_Colm, notekeys="", designD_ncol)
{
	if ( design_Colm > designD_ncol )
	{
		stop(c("The '", notekeys, "' is larger than ",
			   "the total number of columns of designfile, which is ",
			   designD_ncol, ". The reason may be mistakenly setting one or ",
			   "both of these two parameters, or mistakenly setting the 'filestyle'.\n") )
	}
}
######---------------------------------------
##check data type subject by subject
checkDataTypeF <- function(subject_dgl)
{
    nonIntegerV <- unevenV <- missV <- dupV <- uniNumV <- NULL
	#how many NA are permitted to be added in the data
	#here 'maxfoldNA = 3' means fold change between 'timepoints_addedNA' and 
	#original 'timepoints' should be no larger than 3
	#(eg. '0, 2, NA, 6, 8, NA, 12' is permitted, 
	#while '0, NA, NA, NA, 4, NA, NA, NA, 8' is not permitted)
	maxfoldNA <- 3
	for (i in 1:length(subject_dgl))
	{
		tepD <- subject_dgl[[i]]
		tepTime <- sort(tepD$timeT)
		uniTime <- unique(tepTime)
		uniNumV <- c(uniNumV, length(uniTime))
		
		if (length(uniTime) > 1)  
		{
			if ( length(tepTime) != length(uniTime) )
			{	
				dupV <- c(dupV, TRUE)
			}  else  {
				dupV <- c(dupV, FALSE)
			}
			
			uniTimeDiff <- diff(uniTime)
			if ( !all( round(uniTimeDiff) == uniTimeDiff ) )
			{	
				nonIntegerV <- c(nonIntegerV, TRUE)
			}  else  {
				nonIntegerV <- c(nonIntegerV, FALSE)
			}
			
			numDiff <- length(unique(uniTimeDiff))
			minDiff <- min(uniTimeDiff)
			uniFold <- uniTimeDiff/minDiff
			if (numDiff == 1)  {
				unevenV <- c(unevenV, FALSE)
				missV <- c(missV, FALSE)
			}  else if ( !all( round(uniFold) == uniFold ) )  {
				unevenV <- c(unevenV, TRUE)
				missV <- c(missV, FALSE)
			}  else {
				##if too many time points are missing, take it as uneven sampling
				if (length(uniTime)*maxfoldNA <= length(seq(uniTime[1], uniTime[length(uniTime)], by=minDiff)) )  {
					unevenV <- c(unevenV, TRUE)
					missV <- c(missV, FALSE)
				}  else  {
					unevenV <- c(unevenV, FALSE)
					missV <- c(missV, TRUE)
				}
			}
		}  else  {
			nonIntegerV <- c(nonIntegerV, FALSE)
			dupV <- c(dupV, FALSE)
			unevenV <- c(unevenV, FALSE)
			missV <- c(missV, FALSE)
		}
	}
	
	outD <- data.frame("nonInteger"=nonIntegerV, "uneven"=unevenV, 
	                   "miss"=missV, "dup"=dupV, "uniNum"=uniNumV)
	return(outD)
}
######---------------------------------------
##add columns of missing values for 'JTK'
addMissingValueF <- function(dataD, rawtime, duplicate)
{
	unipoints <- unique(rawtime)
	unidiff <- diff(unipoints)
	min_interval <- min(unidiff)
	all_unipoints <- seq(unipoints[1], unipoints[length(unipoints)], by=min_interval)
	if (duplicate) {
		all_dupnum <- rep(-1, length(all_unipoints))
		names(all_dupnum) <- all_unipoints
		input_dupnum <- table(rawtime)
		all_dupnum[names(input_dupnum)] <- input_dupnum
		miss_points <- as.numeric( names(all_dupnum[all_dupnum == -1]) )
		all_dupnum[all_dupnum == -1] <- 1
		all_points <- rep(all_unipoints, as.numeric(all_dupnum) )
		
		if (all(miss_points %in% all_points))  {
			miss_index <- match(miss_points, all_points)
			outD <- as.data.frame(matrix(NA, nrow = nrow(dataD), ncol = length(all_points) ))
			outD[,-(miss_index)] <- dataD
			outname <- paste("misspoint", all_points, sep="_")
			outname[-(miss_index)] <- colnames(dataD)
			outD <- cbind(rownames(dataD), outD)
			dimnames(outD) <- list("r" = rownames(dataD), "c" = c("CycID", outname) )
		}  else  {
			stop(c("There is a unknown bug for 'addMissingValueF()' in 'meta3d', ",
			       "please contact the author.\n"))
		}  
	}  else  {
		all_points <- all_unipoints
		outname <- paste("misspoint", all_points, sep="_")
		names(outname) <- all_points
		outD <- as.data.frame( matrix(NA, nrow = nrow(dataD), ncol = length(all_points) ) )
		colnames(outD) <- all_points
		uniname <- as.character(unipoints)
		outD[,uniname] <- dataD
		outD <- cbind(rownames(dataD),outD)
		outname[uniname] <- colnames(dataD)
		dimnames(outD) <- list("r" = rownames(dataD), "c" = c("CycID", outname))
	}
	return(list("out" = outD, "tim" = all_points))
}
######---------------------------------------	
##use meta2d to analyze time-series datasets subject by subject
##column and row names of "indata" are libraryID and row order in "datafile", respectively
##the "intime" is already sorted in order; "datatype" is a vector containing five elements-
##nonInteger, uneven, miss, dup, uniNum, see the return values of "checkDataTypeF"
runmeta2dF <- function(indata, intime, datatype, rundir, runMethod="JTK", 
                       minper=20, maxper=28, ARSmle="auto", ARSdefaultPer=24)
{
	meta2data <- indata
	meta2time <- intime
	idname <- rownames(indata)
	meta2data <- cbind(idname, meta2data)
	colnames(meta2data) <- c("CycID", colnames(indata) )
	rownames(meta2data) <- idname
	
	datatype <- unlist(datatype[,1:4])
	misstype <- datatype[3]
	duptype <- datatype[4]
	##add 'NA' if there is missing time point and hope to run 'JTK'
	if ( (runMethod == "JTK") & (misstype) )  {
		addNA <- addMissingValueF(dataD=indata, rawtime=intime, duplicate=duptype)
		meta2data <- addNA$out
		meta2time <- addNA$tim
	}

	##write 'meta2data' to a temporary file 
	tempname <- paste("runmeta2dtemp", as.character(runif(1,0,1)), sep="")
	tempfname <- paste(rundir, .Platform$file.sep, tempname, sep="")
	write.table(meta2data, file=tempfname, quote=FALSE, sep="\t", row.names=FALSE)	
	##run selected method with 'meta2d()' function	
	outL <- meta2d(infile=tempfname, outdir=rundir, filestyle="txt", timepoints=meta2time,
                minper=minper, maxper=maxper, cycMethod=runMethod, adjustPhase="predictedPer", 
				combinePvalue="fisher", analysisStrategy="auto", weightedPerPha=FALSE,
                outputFile=FALSE, ARSmle=ARSmle, outIntegration="onlyIntegration", 
                ARSdefaultPer=ARSdefaultPer, outRawData=TRUE, releaseNote=FALSE)
    outD <- outL$meta
	##delete the temporary file and return results
	file.remove(tempfname)
	#file.remove(paste(rundir, .Platform$file.sep, "meta2d_", tempname, sep=""))
	return(outD)
}
######======================================================================================================================================
##integrate multiple P-values by "fisher method"
getCirCasePVA <- function(pvalueM, method="Fisher")
{
    ##the first column in pvalueM is id name, and other columns are p-values
    pvalueID <- pvalueM[,1];
    #pvaM <- matrix(rep(NA,nrow(pvalueM)),ncol=1);
    if ( (method == "Fisher") | (method == "fisher") ) {
		charM <- pvalueM[,-1];
		char2M <- matrix(as.numeric(charM),nrow=nrow(charM),ncol=ncol(charM));
		if ( (!all(as.numeric(charM[,1]) == char2M[,1])) | 
		        (!all(as.numeric(charM[1,]) == char2M[1,])) ) {
			stop("There is a bug in 'getCirOmicsPVA', please contact the author.\n");
		}
		char2M[is.na(char2M) | is.nan(char2M)] <- 1;
		pvalue_meta <- fisher.method(char2M, p.corr="BH", zero.sub = 1e-50);
		pvalue_meta <- as.numeric(pvalue_meta$p.value);
		qvalue <- p.adjust(pvalue_meta, method="BH");
		pvaM <- cbind(pvalue_meta, qvalue);
		rownames(pvaM) <- pvalueID;
		return(pvaM);
	} else {
		stop(c("If the method is correctly set as 'Fisher' or 'fisher', there is ",
		       "unknown bug in 'getCirCasePVA()'. Please contact the author.\n") );
	}
}
######---------------------------------------
##transfer p-value to weight value
getweightM <- function(pvalueM)
{
    weitID <- pvalueM[,1];
    weitM <- apply(pvalueM[,-1], 1, function(z) {
                                 z <- as.numeric(z);
                                 z[is.na(z) | is.nan(z)] <- 1;
                                 ##when Pvalue smaller than 1e-50, set it as 1e-50; different with 'meta2d' (z[!z] <- 1e-300)
                                 z[z < 1e-50] <- 1e-50
                                 z <- -log10(z);
                                 return(z);
                                 });
    weitM <- t(weitM);
    dimnames(weitM)[[1]] <- weitID;
    return(weitM);
}
######---------------------------------------
##average multiple period and phase values
getCirCasePerPha <- function(periodM, phaseM, pvalueM, adjustV, WEIT)
{
    ##the first column of periodM, phaseM and pvalueM is id name (id order) 
    ##other columns are period, phase and p-value information
    perpha_ID <- periodM[,1];
    weitM <- getweightM(pvalueM);
    perphaM <- cbind(periodM[,-1], phaseM[perpha_ID,-1], weitM[perpha_ID,]);
    ADJL <- adjustV;
    caseMF <- function(z, ADJL, WEIT) {
		z <- as.numeric(z);
		zper <- z[1:(length(z)/3)];
		zpha <- z[(length(z)/3+1):(length(z)/3*2)];
		zwei <- z[(length(z)/3*2+1):length(z)];
		if (!WEIT)                                                                 
		{   zwei <- rep(1, length(zwei));    }
		##for JTK, it may output period length as '0', which needs to remove before calculating mean period length
		per_index <- which(!is.na(zper) & !is.nan(zper) & (zper > 0));
		if (length(per_index) > 0)
		{
			zper_mean <- sum(zper[per_index]*zwei[per_index])/sum(zwei[per_index]);
		}  else  {
			zper_mean <- NA;
		}
		pha_index <- which(!is.na(zpha) & !is.nan(zpha));
		##until now, have not found phase was a numeric value, but period is 'NA' or '0'; so zper[pha_index] works
		if ( (length(per_index) > 0) & (length(pha_index) > 0) )
		{
			##the 'circularMean()' function is in 'metaSubF.R'
			zpha_mean <- circularMean(zpha[pha_index], subper=zper[pha_index], 
		                          zweit=zwei[pha_index], meanper=zper_mean, subadj=ADJL);
		}  else  {
			zpha_mean <- NA;
		}	
		return(c(zper_mean,zpha_mean));
	}
    perpha_avgM <- apply(perphaM, 1, caseMF, ADJL=ADJL, WEIT=WEIT);
    perpha_avgM <- t(perpha_avgM);
    per_avg <- perpha_avgM[,1];
    pha_avg <- perpha_avgM[,2];
    ##adjusting the average phase
    pha_avg <- subAdjPha(pha_avg, subper=per_avg, adjustV=ADJL);                                    
    outM <- cbind(per_avg, pha_avg);
    rownames(outM) <- perpha_ID;
    return(outM);
}
######-------------------------------------
##average baseline and relative amplitude
getCirCaseBaseAMP <- function(baseM, ampM, rampM, pvalueM, WEIT)
{
    baseID <- baseM[,1];
    weitM <- getweightM(pvalueM);
    baseweitM <- cbind(baseM[,-1], ampM[baseID, -1], rampM[baseID,-1], weitM[baseID,]);
    caseBF <- function(z, WEIT) {
		z <- as.numeric(z);
		zbase <- z[1:(length(z)/4)];
		zamp <- z[(length(z)/4 + 1): (length(z)/4*2)];
		zramp <- z[(length(z)/4*2 + 1): (length(z)/4*3)];
		zwei <- z[(length(z)/4*3 + 1):length(z)];
		if (!WEIT)                                                         
		{   zwei <- rep(1, length(zwei));    }
		base_index <- which(!is.na(zbase) & !is.nan(zbase));
		base_mean <- sum(zbase[base_index]*zwei[base_index])/sum(zwei[base_index]);
		amp_index <- which(!is.na(zamp) & !is.nan(zamp));
		amp_mean <- sum(zamp[amp_index]*zwei[amp_index]) / sum(zwei[amp_index]);
		ramp_index <- which(!is.na(zramp) & !is.nan(zramp));
		ramp_mean <- sum(zramp[ramp_index]*zwei[ramp_index])/sum(zwei[ramp_index]);
		return(c(base_mean, amp_mean, ramp_mean));
	}
    base_avgM <- apply(baseweitM, 1, caseBF, WEIT=WEIT);
    base_avgM <- t(base_avgM);
    dimnames(base_avgM) <- list("r"=baseID,"c"=c("base", "amp", "ramp"));
    return(base_avgM);
}
######---------------------------------------
##integrate analysis results from multiple subjects
##and output integrated P-value, period, phase, baseline, AMP and rAMP values
integrateLSOUT <- function(lsoutL, adjustV, weightS)
{
    lsLen <- length(lsoutL);
    WEIT <- weightS;
    ##check whether lsoutL is empty
    if (lsLen)
    {
        initD <- lsoutL[[1]];
        ##the column order of initD: 'CycID', 'Pvalue', 'BH.Q', 'Period', 'Phase',  
        ##and 'Amplitude', 'meta2d_Base', 'meta2d_AMP', 'meta2d_rAMP', with or without raw values
        initD <- initD[,1:9];
        generalID <- as.character(initD[,"CycID"]);
        header <- c("CycID", "meta3d_Pvalue", "meta3d_BH.Q", "meta3d_Period", 
                    "meta3d_Phase", "meta3d_Base", "meta3d_AMP", "meta3d_rAMP");
        ##do not use calculated amplitude by each method
        outD <- initD[,-6];                    
        dimnames(outD) <- list("r"=generalID, "c"=header);
		
        if (lsLen >= 2)
        {
            ##use 'Pvalue', 'Period', 'Phase', 'Base', 'AMP', and 'rAMP'
            ##do not use 'ID'(col 1), 'BH.Q'(col 3) and 'Amplitude'(col 6) column
            integrateD <- initD[, -c(1,3,6)];
            dimnames(integrateD)[[1]] <- generalID;
			
            flag <- NULL
            for (i in 2:lsLen)
            {
                processD <- lsoutL[[i]];
                processD <- processD[,1:9];
                ##check whether have same ID order
                if (all(generalID == as.character(processD[,"CycID"])) )
                {
                    integrateD <- cbind(integrateD, processD[, -c(1,3,6)]);
                }  else  {
                    flag <- c(flag, i);
                }
            }
            if (length(flag) > 0)
            {   cat("Not all subjects have the same row names. Please check the input 'datafile'.\n");    }
            
            ##integrating p-value, period, phase, baseline, amplitude and relative amplitude
            if (ncol(integrateD) > 6)
            {
                ##divide integrateD into multiple sub-matrix
                subjectNum <- ncol(integrateD)/6;
                pvalueM <- periodM <- phaseM <- baseM <- ampM <- rampM <- generalID;
                for (k in 1:subjectNum)
                {
                    pvalueM <- cbind(pvalueM, integrateD[,(k-1)*6+1]);
                    periodM <- cbind(periodM, integrateD[,(k-1)*6+2]);
                    phaseM <- cbind(phaseM, integrateD[,(k-1)*6+3]);
                    baseM <- cbind(baseM, integrateD[,(k-1)*6+4]);
                    ampM <- cbind(ampM, integrateD[,(k-1)*6+5]);
                    rampM <- cbind(rampM, integrateD[,(k-1)*6+6]);
                }
                dimnames(pvalueM)[[1]] <- dimnames(periodM)[[1]] <- dimnames(phaseM)[[1]] <- 
                dimnames(baseM)[[1]] <- dimnames(ampM)[[1]] <- dimnames(rampM)[[1]] <- generalID;
                
                ##call associated integration function
                cirhomoPvaM <- getCirCasePVA(pvalueM, method="Fisher");                                  ##'pvalue_meta', 'qvalue'
                cirhomoPerPhaM <- getCirCasePerPha(periodM, phaseM, pvalueM, adjustV, WEIT);             ##'per_avg', 'pha_avg'
                cirhomoBaseAmpM <- getCirCaseBaseAMP(baseM, ampM, rampM, pvalueM, WEIT);                 ##'base', 'AMP', 'rAMP'
                
                ##replace associated columns with integrated values
                ##an internal function
                transferValueF <- function(dataD, cname, dataM, rname)  {
                    if (nrow(dataD) == nrow(dataM))  {
                        for (k in 1:ncol(dataM))
                        {  dataD[, cname[k]] <- dataM[rname,k];  }
                        return(dataD);
                    }  else  {
                         stop("There is bug in 'integrateLSOUT()', please contact the author.\n");
                    }
                }
                outD <- transferValueF(dataD=outD, cname=c("meta3d_Pvalue","meta3d_BH.Q"), dataM=cirhomoPvaM, rname=generalID);
                outD <- transferValueF(dataD=outD, cname=c("meta3d_Period","meta3d_Phase"), dataM=cirhomoPerPhaM, rname=generalID);
                outD <- transferValueF(dataD=outD, cname=c("meta3d_Base","meta3d_AMP","meta3d_rAMP"), 
                                       dataM=cirhomoBaseAmpM, rname=generalID);
           }  else  {
              ##adjust phase; column order of 'integrateD'-
              ##the 'Pvalue', 'Period', 'Phase', 'Base', 'AMP', and 'rAMP'
              cirhomoPha <- subAdjPha(subpha=integrateD[,3], subper=integrateD[,2], adjustV);
              outD[,"meta3d_Phase"] <- cirhomoPha;
           }
        }  else if  (lsLen == 1)  {
            ##adjust phase; column order of 'outD'-'CycID', 'meta3d_Pvalue', 'meta3d_BH.Q', 
            ##the 'meta3d_Period', 'meta3d_Phase', 'meta3d_Base', 'meta3d_AMP', 'meta3d_rAMP'
            cirhomoPha <- subAdjPha(subpha=outD[,"meta3d_Phase"], subper=outD[,"meta3d_Period"], adjustV);
            outD[,"meta3d_Phase"] <- cirhomoPha;
        }  else  {
            stop("There is unknown bug associated with 'integrateLSOUT()'. Please contact the author.\n");
        }
        return(outD);
    }  else  {
        cat(c("There is no integration result for one or more groups, please check ",
		      "whether each parameter is correctly setted according to instruction.\n"));
        return(NA);
    }
}
