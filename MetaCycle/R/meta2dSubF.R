### This file contains functions used in 'meta2dMainF.R'
### Author: Gang Wu(wggucas@gmail.com)
### Lab: John Hogenesch's lab (http://hogeneschlab.org/)
###======================================================================================================================================
##check timepoints
checkTimeF <- function(timepoints)
{
	NUMT <- TRUE;
	if (is.character(timepoints))
    {
        if ( (timepoints == "Line1") | (timepoints == "line1") )
        {   NUMT <- FALSE; }
        else
        {
            stop(c("If the time points values are in the first line of ",
               "infile, please set 'timepoints' as 'Line1' or 'line1'.\n") );
        }
    }  else if ( is.numeric(timepoints) )  {
         if ( any(is.na(timepoints)) | any(is.nan(timepoints)) )
         {  stop("The 'timepoints' should not contain 'NA' or 'NaN'.\n");    }
         if (min(timepoints) < 0 )
         {
            stop("The 'timepoints' should be a non-negative numeric vector.\n");
         }
    }  else {
        stop(c("The 'timepoints' should be set as a non-negative numeric vector ",
           "value. If the first line of 'infile' contains timepoints ",
           "information, 'timepoints' could be set as 'Line1' or 'line1'.\n") );
    }
	return(NUMT);
}
##-------------------------------------
##check 'ARSdefaultPer'
checkARSdefaulterPerF <- function(ARSdefaultPer, minper, maxper)
{
	if (length(ARSdefaultPer) == 1)
	{
		if (is.numeric(ARSdefaultPer))
		{
			ARS_PER <-  ARSdefaultPer;
			if ( (ARS_PER < minper) | (ARS_PER > maxper) )
			{
				stop(c("Please set 'ARSdefaultPer' as the expected ",
				  "period length, but it should be no smaller than ",
				  "'minper' and no larger than 'maxper'. \n") );
			}
		}  else  {
			stop(c("Please give one numeric value to 'ARSdefaultPer'. ",
				"It is suggested to set as the expected period length.\n") );
		}
	} else {
		stop(c("Please give one numeric value to 'ARSdefaultPer'. ",
			"It is suggested to set as the expected period length.\n") );
	}
	return(ARS_PER);
}
##-------------------------------------
##adjust phase values from ARS with defined period length, and also 
##return unique phase and period value for each profile  
adjPhaARS <- function(arsM,adjustV)
{
    ##the 'arsM' is a matrix containing all the output values from ARS, 
    ##adjustV is the defined period length used for adjusting phase
    if (adjustV == 0) {
        ##if the 'adjustV' is set '0', output the original values from ARS 
        return(arsM[,c("period","phase","amplitude")]);                                            
    } else {
        arsid <- dimnames(arsM)[[1]];
        arsnum <- as.numeric(arsM[,"period_number"]);
        names(arsnum) <- arsid;
        arsnum[is.na(arsnum) | is.nan(arsnum)] <- 0;
        uni_id <- names(arsnum[arsnum <= 1]);
        multi_id <- names(arsnum[arsnum > 1]);
        ##select the period and phase information corresponding to 
        ##the largest amplitude if multiple periods exist
        if (length(multi_id) > 0) {
            catM <- arsM[multi_id,c("period","phase","amplitude")];
            if (length(multi_id) == 1)
            { catM <- matrix(catM, nrow=1, ncol=length(catM));  }
            outM <- apply(catM,1,function(z){
                            zper <- unlist(strsplit(z[1],","));
                            names(zper) <- 1:length(zper);
                            zpha <- unlist(strsplit(z[2],","));
                            names(zpha) <- 1:length(zpha);
                            zamp <- unlist(strsplit(z[3],","));
                            zamp <- as.numeric(zamp);
                            names(zamp) <- 1:length(zamp);
                            zamp<- sort(zamp,decreasing=TRUE);
                            maxid<- names(zamp[1]);
                            return(c(zper[maxid],zpha[maxid],zamp[maxid]));
                        });
            ##each column corresponding to one profile in the results from 'apply(,1,)' 
            outM <- t(outM);                                                                        
            outM <- rbind(outM,arsM[uni_id,c("period","phase","amplitude")]);
            dimnames(outM) <- list("r"=c(multi_id,uni_id),"c"=c("period","phase","amplitude"));
            outM <- outM[arsid,];
        } else {
            outM <- arsM[arsid,c("period","phase","amplitude")];
        }
        ##adjust phase value with defined period length
        adj_outM <- outM;
        pha <- as.numeric(adj_outM[,"phase"]);
        per <- as.numeric(adj_outM[,"period"]);
        adj_pha <- subAdjPha(pha, subper=per, adjustV);
        adj_outM <- cbind(adj_outM[,"period"],adj_pha,adj_outM[,"amplitude"]);
        dimnames(adj_outM) <- list("r"=arsid,"c"=c("period","phase","amplitude"));
        return(adj_outM);
    }
}
##-------------------------------------
##adjust output phase values from JTK with defined period length
adjPhaJTK <- function(jtkM, adjustV, START_TIME)
{
    ##the 'jtkM' is a matrix containing all the output values from JTK
    jtkid <- dimnames(jtkM)[[1]];
    pha <- as.numeric(jtkM[,"LAG"]);
    if (adjustV == 0) {
        adjpha <- pha;
    }  else  {
        pha <- pha + START_TIME;
        per <- as.numeric(jtkM[,"PER"]);
        adjpha <- subAdjPha(pha, subper=per, adjustV);
    }
    names(adjpha) <- jtkid;
    return(adjpha);
}
##-------------------------------------
##adjust output phase value from LS with defined period length
adjPhaLS <- function(lsM,adjustV)
{
    ##the 'lsM' is a matrix containing all the output values from LS
    lsid <- dimnames(lsM)[[1]];
    pha <- as.numeric(lsM[,"PhaseShift"]);
    if (adjustV == 0) {
        adjpha <- pha;
    }  else  {
        per <- as.numeric(lsM[,"Period"]);
        adjpha <- subAdjPha(pha, subper=per, adjustV);
    }
    names(adjpha) <- lsid;
    return(adjpha);
}
##-------------------------------------
##output separate analysis results
outSigResultF <- function(outM, SIG_row, outRawData, EXPM, ID_ORDER, methodName, 
                          outfile_name,   features,   tabletypeL,   outputFile)
{
	ONE_OUTM <- outM
	if (length(ONE_OUTM) > 1) {
		if (outputFile)
		{
			if (!SIG_row)
			{
				if (!outRawData)
				{
					write.table(ONE_OUTM,file=outfile_name[methodName],row.names=FALSE,
								sep=tabletypeL$sep,quote=tabletypeL$quote,dec=tabletypeL$dec);
				}  else  {
					write.table(cbind(ONE_OUTM[ID_ORDER,],EXPM[ID_ORDER,]),file=outfile_name[methodName],
								row.names=FALSE,sep=tabletypeL$sep,quote=tabletypeL$quote,dec=tabletypeL$dec);
				}
			}  else  {
				if (!outRawData)
				{
					sigoutM <- data.frame(ONE_OUTM, stringsAsFactors=FALSE);
				}  else  {
					sigoutM <- data.frame(cbind(ONE_OUTM[ID_ORDER,], EXPM[ID_ORDER,]), stringsAsFactors=FALSE);
				}
				##keep 'fdr' value and 'pvalue' same if 'SIG_row=TRUE'
				sigoutM[,features[1]] <- sigoutM[,features[2]];
				write.table(sigoutM[1,],file=outfile_name[methodName],row.names=FALSE,
							sep=tabletypeL$sep,quote=tabletypeL$quote,dec=tabletypeL$dec);
			}
		}
	} else {
		stop(c("Warning: no output result from ", methodName, ". ",
		     "If 'analysisStrategy' was set as 'selfUSE', please check whether ",
			 "it is suitable for analyzing the input file. If 'analysisStrategy' ",
		     "was set as 'auto', please contact the author.\n") );
	}
}
###======================================================================================================================================
##get the integrated P-value from P-values calculated by multiple methods
getCirOmicsPVA <- function(pvalueM, method="bonferroni")
{
    ##the first column in pvalueM is id name, and other columns are p-values
    pvalueID <- pvalueM[,1];
    pvaM <- matrix(rep(NA,nrow(pvalueM)),ncol=1);
    if ( (method == "bonferroni") | (method == "Bonferroni") ) {
		pvalue_meta <- NA;
		pvalue_meta <- apply(pvalueM[,-1],1,function(z) {
              z <- as.numeric(z);
              z[is.na(z) | is.nan(z)] <- 1;
              z_meta <- p.adjust(z, method="bonferroni");
              z <- min(z_meta);
              return (z);
						   } )
		qvalue <- p.adjust(pvalue_meta, method="BH");
		pvaM <- cbind(pvalue_meta, qvalue);
	} else if ( (method == "Fisher") | (method == "fisher") ) {
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
	} else {
		stop(c("If the method is correctly set as 'bonferroni', 'Bonferroni', 'fisher' or ",
		     "'Fisher', there is unknown bug in combining P-values. Please contact the author.\n") );
	}
	rownames(pvaM)<- pvalueID;
	return (pvaM);
}
##-------------------------------------
##get the average period and phase
getCirOmicsPerPha <- function(periodM, phaseM, pvalueM, adjustV, weightedPerPha)
{
    ##the first column of periodM, phaseM and pvalueM is id name
    perpha_ID <- periodM[,1];
    ADJ <- adjustV;
    ##prepare weight matrix
    rankM <- apply(pvalueM[perpha_ID,-1], 2, function(z) {
                                 z <- as.numeric(z);
                                 z[is.na(z) | is.nan(z)] <- 1;
                                 z[!z] <- 1e-300;
                                 z <- -log10(z);
                                 return (rank(z));
                                 });
    perphaM <- cbind(periodM[,-1], phaseM[perpha_ID,-1], rankM);
    ##an internal function for calculating average period length and phase value
    zmeanF <- function(z, ADJ, weightedPerPha)  {
		z <- as.numeric(z);
		zper <- z[1:(length(z)/3)];
		zpha <- z[(length(z)/3+1):(length(z)/3*2)];
		zwei <- z[(length(z)/3*2+1):length(z)];
		if (!weightedPerPha)
		{   zwei <- rep(1, length(zwei));    }
		##for JTK, it may output period length as '0', which needs to remove before calculating mean period length
		per_index <- which(!is.na(zper) & !is.nan(zper) & (zper > 0) );
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
			zpha_mean <- circularMean(zpha[pha_index], subper=zper[pha_index], 
		                          zweit=zwei[pha_index], meanper=zper_mean, subadj=ADJ);
		}  else  {
			zpha_mean <- NA;
		}
		return(c(zper_mean, zpha_mean));
	}
    perpha_avgM <- apply(perphaM, 1, zmeanF, ADJ=ADJ, weightedPerPha=weightedPerPha);
    perpha_avgM <- t(perpha_avgM);
    per_avg <- perpha_avgM[,1];
    pha_avg <- perpha_avgM[,2];
    ##adjust the average phase
    pha_avg <- subAdjPha(pha_avg, subper=per_avg, ADJ);  
    outM <- cbind(per_avg, pha_avg);
    rownames(outM) <- perpha_ID;
    return(outM);
}
##-------------------------------------
##calculate baseline and amplitude using ordinary least square method(OLS)
##set a general expression function as: Y = Baseline + Trend*(tim - mean(tim)) + Amplitude*cos(2*pi/Period*(tim - Phase))
##Period and Phase are calculated by 'getCirOmicsPerPha', using OLS to estimate three unknown parameters,
##Baseline(Base), Trend and Amplitude(AMP). The relative Amplitude(rAMP) is defined as Amplitude/Baseline.
getCirOmicsAMP <- function(exprperphaM, AMPTIM)
{
    ##an internal function for calculating amplitude
    getAMP <- function(z, AMPTIM) {
		z <- as.numeric(z);
		pha <- z[length(z)];
		per <- z[length(z) - 1];
		expr <- z[1:(length(z) - 2)];
		##removing time points corresponding to missing values
		timv <- tim <- AMPTIM[!is.na(expr) & !is.nan(expr)];
		exprv <- expr <- expr[!is.na(expr) & !is.nan(expr)];
		trendt <- tim - mean(tim[!is.na(tim) & !is.nan(tim)]);
		if ( (is.na(per)) | (is.na(pha)) | is.nan(per) | is.nan(pha) )
		{
			basev <- median(exprv);
			out <- c(basev, NA, NA);
		}  else  {
			cost <- cos(2*pi/per*(tim - pha));
			fit <- lm(expr~trendt + cost);
			fitcoef <- fit$coefficients;
			basev <- fitcoef[1];
			ampv <- fitcoef[3];
			if ( (basev < min(exprv)) | (basev > max(exprv)) ) {
				basev <- median(exprv);
				##calculate the mean value if with replicates; if '3.1, 2.9' are two replicate values in one time point, 
				##the 'detrend_linear()' will give different results when changing the order to '2.9, 3.1'
				if (length(exprv) > length(unique(timv)) )                        
				{
					timv_uni <- sort(unique(timv));
					timv <- factor(timv, levels=timv_uni);
					exprv_mean <- tapply(exprv, timv, mean);
					exprv_mean <- as.numeric(exprv_mean);
					timv <- timv_uni;
					exprv <- exprv_mean;
				}
				##the 'detrend_linear()' is in 'ARS.R', detrending the obvious trend in the profile
				exprv_dt <- detrend_linear(exprv);                                     
				costv <- cos(2*pi/per*(timv - pha));
				fitv <- lm(exprv_dt~costv);
				fitcoefv <- fitv$coefficients;
				ampv <- fitcoefv[2];
			} 
			if (abs(basev) >= 1) { 
				out <- c( basev, abs(ampv), abs(ampv/basev) );
			} else {
				out <- c( basev, abs(ampv), abs(ampv) );
			}
		}
		names(out) <- c("Baseline", "Amplitude", "RelativeAmplitude");
		return(out);
	}
	##each column corresponding to one profile in the results from 'apply(,1,)' 
	outM <- apply(exprperphaM, 1, getAMP, AMPTIM=AMPTIM);
	outM <- t(outM);                                                                               
  return(outM);
}
##-------------------------------------
##get integration header
getSigIntegHeaderF <- function(methodName, INTEGRATION, integration_header)
{
	sub_header <- c("pvalue","BH.Q","period","adjphase","amplitude");
	if (INTEGRATION == "onlyIntegration")
	{
		integration_header <- c(integration_header,paste(methodName, sub_header, sep="_"));
	}
	integration_header <- c(integration_header, "meta2d_Base", "meta2d_AMP", "meta2d_rAMP");
	return(integration_header);	
}
##-------------------------------------
##get integration result when only one method is used
getSigIntegResultF <- function(sigoutL, outfile_tag, EXPM, ID_ORDER, ADPHA, AMPTIM, 
                               START_TIME, INTEGRATION, integration_outM, integration_header)
{
	ARS_OUTM <- sigoutL$ARS;
	JTK_OUTM <- sigoutL$JTK;
	LS_OUTM <- sigoutL$LS;
	if ( (outfile_tag["ARS"]) & (length(ARS_OUTM) > 1) )
	{
		ars_adjM <- adjPhaARS(ARS_OUTM[ID_ORDER,],adjustV=ADPHA);
		arsEXPMPERPHA <- cbind(EXPM[ID_ORDER,],ars_adjM[ID_ORDER,1:2]);
		arsampM <- getCirOmicsAMP(exprperphaM=arsEXPMPERPHA, AMPTIM);
		rownames(arsampM) <- ID_ORDER;
		integration_outM <- cbind(integration_outM, ARS_OUTM[ID_ORDER,c("pvalue","fdr_BH")],
								  ars_adjM[ID_ORDER,], arsampM[ID_ORDER,]);
		integration_header <- getSigIntegHeaderF(methodName="ARS", INTEGRATION, integration_header);
	}
	if ( (outfile_tag["JTK"]) & (length(JTK_OUTM) > 1) )
	{
		jtk_adjpha <- adjPhaJTK(JTK_OUTM[ID_ORDER,], adjustV=ADPHA, START_TIME);
		jtkEXPMPERPHA <- cbind(EXPM[ID_ORDER,], JTK_OUTM[ID_ORDER,"PER"], jtk_adjpha[ID_ORDER]);
		jtkampM <- getCirOmicsAMP(exprperphaM=jtkEXPMPERPHA, AMPTIM);
		rownames(jtkampM) <- ID_ORDER;
		integration_outM <- cbind(integration_outM, JTK_OUTM[ID_ORDER,c("ADJ.P","BH.Q","PER")],
								  jtk_adjpha[ID_ORDER], JTK_OUTM[ID_ORDER,"AMP"], jtkampM[ID_ORDER,]);
		integration_header <- getSigIntegHeaderF(methodName="JTK", INTEGRATION, integration_header);
	}
	if ( (outfile_tag["LS"]) & (length(LS_OUTM) > 1) )
	{
		ls_adjpha <- adjPhaLS(LS_OUTM[ID_ORDER,],adjustV=ADPHA);
		lsEXPMPERPHA <- cbind(EXPM[ID_ORDER,], LS_OUTM[ID_ORDER, "Period"], ls_adjpha[ID_ORDER]);
		lsampM <- getCirOmicsAMP(exprperphaM=lsEXPMPERPHA, AMPTIM);
		rownames(lsampM) <- ID_ORDER;
		integration_outM <- cbind(integration_outM, LS_OUTM[ID_ORDER,c("p","BH.Q","Period")],
								  ls_adjpha[ID_ORDER], LS_OUTM[ID_ORDER,"PhaseShiftHeight"], lsampM[ID_ORDER,]);
		integration_header <- getSigIntegHeaderF(methodName="LS", INTEGRATION, integration_header);
	}
	sigIntegL <- list("integ"=integration_outM, "header"=integration_header);
	return(sigIntegL);
}
##-------------------------------------
##output integrated analysis results
outIntegResultF <- function(integration_outM, EXPM, ID_ORDER, integration_header, 
                            SIG_row, outRawData, integ_outname, tabletypeL)
{
	if (!SIG_row)
	{
		if (!outRawData)
		{
			write.table(integration_outM, file=integ_outname, row.names=FALSE,
			            sep=tabletypeL$sep, quote=tabletypeL$quote, dec=tabletypeL$dec);
		}  else  {
			write.table(cbind(integration_outM[ID_ORDER,], EXPM[ID_ORDER,]),file=integ_outname,
			                 row.names=FALSE, sep=tabletypeL$sep, quote=tabletypeL$quote, dec=tabletypeL$dec);
		}
	}  else  {
		if (!outRawData)
		{
			integration_outM <- data.frame(integration_outM,stringsAsFactors=FALSE);
		}  else  {
			integration_outM <- data.frame(cbind(integration_outM[ID_ORDER,], EXPM[ID_ORDER,]), stringsAsFactors=FALSE);
		}
		pva_index <- grep("_pvalue", integration_header);
		qva_index <- grep("_BH.Q", integration_header);
		integration_outM[,qva_index] <- integration_outM[,pva_index];
		write.table(integration_outM[1,], file=integ_outname, row.names=FALSE,
		            sep=tabletypeL$sep, quote=tabletypeL$quote, dec=tabletypeL$dec);
	}
}
##-------------------------------------
##return analysis results if 'outputFile=FALSE'
returnResultF <- function(outLIST, integration_outM, EXPM, ID_ORDER, outfile_tag, SIG_row, outRawData)
{
	##an internal function
	getNumericF <- function(inputM, expM, id_order, outRaw, methodName, sigRow, initcol=1)  {
		header <- dimnames(inputM)[[2]];
		if (outRaw)
		{  
			inputM <- cbind(inputM[id_order,], expM[id_order,]);
			header <- c(header, dimnames(expM)[[2]]);
		}
		inputD <- as.data.frame(inputM, stringsAsFactors=FALSE);
		dimnames(inputD)[[1]] <- id_order;
		dimnames(inputD)[[2]] <- header;
		if (ncol(inputD) >= initcol)
		{
			for (i in initcol:ncol(inputD)) 
			{
				inputD[,i] <- as.numeric(inputD[,i]);
			}
		}  else  {
			stop("There is a bug in 'returnResultF()', please contact the author.\n");
		}
		if (sigRow)  {
			##keep BH.Q and pvalue as same
			if (methodName == "ARS")  {
				inputD[,"fdr_BH"] <- inputD[,"pvalue"];
			}  else if (methodName == "JTK")  {
				inputD[,"BH.Q"] <- inputD[,"ADJ.P"];
			}  else if (methodName == "LS")  {
				inputD[,"BH.Q"] <- inputD[,"p"];
			}  else if (methodName == "meta")  {
				pva_index <- grep("_pvalue", colnames(inputD));
				qva_index <- grep("_BH.Q", colnames(inputD));
				inputD[,qva_index] <- inputD[,pva_index];
			}  else {
				stop("There is a bug in 'returnResultF()' function. Please contact the author.\n");
			}
			return(inputD[1,]);
		}  else  {
			return(inputD);
		}
	}
	##the initial return list
	meta2dL <- list("ARS"=NULL, "JTK"=NULL, "LS"=NULL, "meta"=NULL);
	ARS_OUTM <- outLIST$ARS;
	JTK_OUTM <- outLIST$JTK;
	LS_OUTM <- outLIST$LS;
	##get the return result one by one
	if ( (outfile_tag["ARS"]) & (length(ARS_OUTM) > 1) )
	{
		arsD <- getNumericF(inputM=ARS_OUTM, expM=EXPM, id_order=ID_ORDER, 
		           outRaw=outRawData, methodName="ARS", sigRow=SIG_row, initcol=8);
		meta2dL$ARS <- arsD;
	}
	if ( (outfile_tag["JTK"]) & (length(JTK_OUTM) > 1) )
	{
		jtkD <- getNumericF(inputM=JTK_OUTM, expM=EXPM, id_order=ID_ORDER, 
		           outRaw=outRawData, methodName="JTK", sigRow=SIG_row, initcol=2);
		meta2dL$JTK <- jtkD;
	}
	if ( (outfile_tag["LS"]) & (length(LS_OUTM) > 1) )
	{
		lsD <- getNumericF(inputM=LS_OUTM, expM=EXPM, id_order=ID_ORDER, 
		           outRaw=outRawData, methodName="LS", sigRow=SIG_row, initcol=2);
		meta2dL$LS <- lsD;
	}
	if (length(integration_outM) > 1)
	{
		metaD <- getNumericF(inputM=integration_outM, expM=EXPM, id_order=ID_ORDER, 
					   outRaw=outRawData, methodName="meta", sigRow=SIG_row, initcol=2);
		meta2dL$meta <- metaD;
	}
	return(meta2dL);
}
