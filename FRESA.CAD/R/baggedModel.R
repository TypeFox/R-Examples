baggedModel <-
function(modelFormulas,data,type=c("LM","LOGIT","COX"),Outcome=NULL,timeOutcome=NULL,pvalue=0.05,backElimination=FALSE,frequencyThreshold=0.05,removeOutliers=4.0)
{
	type <- match.arg(type)

	
	loops <- length(modelFormulas);
	theoutcome <- data[,Outcome];
	binoutcome <- (length(table(theoutcome))==2) && (min(theoutcome)==0);
	predtype="linear";
	if (binoutcome) predtype="prob";
	if ( (type=="LM") && (binoutcome==TRUE) )
	{
		data[,Outcome] = 2*theoutcome-1.0;
	}
	features <- vector();
	observations <- nrow(data);
	
#	cat(removeOutliers,"\n");

	for (n in 1:loops)
	{
		frm <- gsub("+ ","+",modelFormulas[n],fixed = TRUE);
		frm <- gsub(" +","+",frm,fixed = TRUE);
		frm <- gsub("  ","",frm,fixed = TRUE);
#		print(frm);
		if (gregexpr(pattern ='~',frm)[1]>0)
		{
			feat <- unlist(strsplit(frm,"[~]"));
			feat <- unlist(strsplit(feat[2],"[+]"));
		}
		else
		{
			feat <- unlist(strsplit(frm,"[+]"));
		}
		for (i in 1:length(feat))
		{
			if ((feat[i]!="1")&&(feat[i]!=" ")&&(feat[i]!="  ")&&(feat[i]!="")) features <- append(features,feat[i]);
		}
	}
	VarFrequencyTable <- table(features);
	VarFrequencyTable <- VarFrequencyTable[order(-VarFrequencyTable)]
	vnames <- rownames(VarFrequencyTable);

	
	nsize <- nrow(data)
	
	baseForm = Outcome;
#For Cox  models 
	if (type == "COX")
	{
	  baseForm = paste("Surv(",timeOutcome);
	  baseForm = paste(baseForm,paste(",",paste(Outcome,")")));
	}

	lastTopVariable = length(VarFrequencyTable);
	if (lastTopVariable >= nsize/3) lastTopVariable <- nsize/3;
	
	frma <- paste(baseForm,"~ ");
	firstvar = ".";
	enterlist <- vector();
	for ( i in 1:length(VarFrequencyTable))
	{
		if ((vnames[i] != " ") && (vnames[i] != ""))
		{
			enterlist <- append(enterlist,vnames[i]);
			if ((i<=lastTopVariable)&&(VarFrequencyTable[i] > frequencyThreshold*loops))  # Only features with a given frequency
			{
				frma <- paste(frma,"+",vnames[i]);
				
				if (firstvar == ".") 
				{
					firstvar = vnames[i];
				}
				
			}
		}
	}
	model <- modelFitting(formula(frma),data,type)
	
	if (binoutcome)
	{
		shortVarList <- enterlist;
		varlist <- vector();
		for (i in 1:length(shortVarList))
		{
			varlist <- append(varlist,str_replace_all(unlist(strsplit(
							str_replace_all(
								str_replace_all(
									str_replace_all(
										str_replace_all(shortVarList[i],"I\\("," ")
									,"\\("," ")
								,">","\\*")
							,"<","\\*")
					,"\\*"))[1]," ",""))
		}
		shortVarList <- c(Outcome,as.vector(rownames(table(varlist))));
#		print(shortVarList)
		casesample = subset(data[,shortVarList],get(Outcome)  == 1);
		controlsample = subset(data[,shortVarList],get(Outcome) == 0);
		nrowcase <- nrow(casesample);
		nrowcontrol <- nrow(controlsample);
		minsample <- min(nrowcase,nrowcontrol);
		if (nrowcase==nrowcontrol) binoutcome=FALSE;
	}
	
	
#	print(summary(model));
	msize <- length(model$coefficients)
	wts <- 0;
	zvalues <- model$coefficients;
	zwts <- model$coefficients;

	for (i in 1:length(model$coefficients))
	{
		model$coefficients[i] <- 0;
		zvalues[i] = 0;
		zwts[i] = 0;
	}
	
	if (type=="COX") 
	{
		zcol=4;
	}
	else
	{
		zcol=3;
	}
	
#	cat("\n");
	avgsize = 0;
	varoutcome <- var(theoutcome);
	totresidual <- rep(0,length(theoutcome));
	contresidual <- totresidual;
	for (n in 1:loops)
	{
	
		if ((n %% 10) == 0) cat(".");
		if (gregexpr(pattern ='~',modelFormulas[n])[1]>0)
		{
			ftmp <- formula(modelFormulas[n]);
		}
		else
		{
			ftmp <- formula(paste(baseForm,"~",modelFormulas[n]));
		}
		if (binoutcome) 
		{
			BDataSet <- rbind(casesample[sample(1:nrowcase,minsample,replace=FALSE),],controlsample[sample(1:nrowcontrol,minsample,replace=FALSE),])
			out <- modelFitting(ftmp,data,type);
		}
		else
		{
			out <- modelFitting(ftmp,data,type);
		}
		if (backElimination)
		{
			if (type=="LM")
			{
				out <- backVarElimination_Res(out,pvalue=pvalue,Outcome=Outcome,data=data,type=type,testType="Ftest")$back.model;				
			}
			else
			{
				out <- backVarElimination_Bin(out,pvalue=pvalue,Outcome=Outcome,data=data,type=type,selectionType="zIDI")$back.model;
			}
		}
		if (!inherits(out, "try-error")) 
		{
			class(model) <- class(out)		
			osize <- length(out$coefficients)
			avgsize = avgsize+osize;
			
			if (osize>0)
			{				
				outcoef  <- summary(out)$coefficients;
				curprediction <- predictForFresa(out,data,predtype)
				residual <- abs(curprediction-theoutcome);
				varresidual <- mean(residual^2);
				w = (varoutcome-varresidual)/varoutcome;
				if (predtype=='prob') 
				{
					w = 2.0*(mean((curprediction>0.5) == theoutcome) - 0.5); # we will use accuracy for binary outcomes
				}
				if (!is.na(w))
				{
					if (w>0) 
					{
						w = w-(1.0-w)*(osize-1.0)/(nsize-osize); # adjusting for model size.
						if (w<0.1) w= 0.1; # we will keep the minumum weight to 0.1
						totresidual <- totresidual+w*residual;
						contresidual <- contresidual + w;
						wts = wts + w;
						for (i in 1:msize)
						{
							for (j in 1:osize)
							{
								if (names(model$coefficients)[i] == names(out$coefficients)[j])
								{
									model$coefficients[i] <- model$coefficients[i] + w*out$coefficients[j]; 
									zvalues[i] <- zvalues[i] + outcoef[j,zcol];
									zwts[i] = zwts[i] + 1;
								}
							}
						}
					}
				}
			}
		}
	}
	avgsize = avgsize/loops;
	cat(". End \n");
	if( wts>0)
	{
			model$coefficients <- model$coefficients/wts;
	}
	for (i in 1:msize)
	{
		if( zwts[i]>0)
		{
			zvalues[i] <- zvalues[i]/zwts[i];
		}
	}
#	print(summary(model));
	baggedPredict <- predictForFresa(model,data,"linear");
	ndata <- as.data.frame(cbind(data[,Outcome],baggedPredict));
	if (predtype == 'linear') # for linear predictior we recalibrate the coefficients for minimum RMSE
	{
		out <- modelFitting(formula("V1~baggedPredict"),ndata,"LM");
		model$coefficients <- model$coefficients*out$coefficients[2];
		model$coefficients[1] <- out$coefficients[1]+model$coefficients[1];
	}
	else
	{
		if (type == 'LOGIT') # for logit predict we set the offset to maximum AUC
		{
#			out <- modelFitting(formula("V1~baggedPredict"),ndata,type);
#			model$coefficients <- model$coefficients*out$coefficients[2];
#			baggedPredict <- predictForFresa(model,data,"linear");
			ndata <- as.data.frame(cbind(data[,Outcome],baggedPredict));
			ndata <- ndata[order(baggedPredict),];
			ncases <- sum(data[,Outcome]);
			ncontrol <- nsize-ncases;
			mauc <- 0;
			sen <- 1;
			spe <- 0;
			threshold <- 1;
			for (i in 1:nsize)
			{	
				event <- 1.0*(ndata[i,1]==1);
				nevent <- 1.0*(event==0);
				sen <- sen - event/ncases;
				spe <- spe + nevent/ncontrol;
				auc <- (sen+spe)/2;
				if (auc>mauc)
				{
					mauc <- auc;
					threshold <- i;
				}
			}
			indx <- threshold;
			preindx <- threshold-1;
			postindx <- threshold+1;
			if (preindx<1) preindx=1;
			if (postindx>nsize) postindx=nsize;
			model$coefficients[1] <- model$coefficients[1]-(ndata[preindx,2]+2.0*ndata[indx,2]+ndata[postindx,2])/4;
#			baggedPredict <- predictForFresa(model,data,"linear");
#			ndata <- as.data.frame(cbind(data[,Outcome],baggedPredict));
#			ndata <- ndata[order(baggedPredict),];
#			print(ndata);
		}
	}


#	print(summary(model));
	reducededset <- data;
	baggedPredict <- predictForFresa(model,data,predtype);
	residual <- baggedPredict-theoutcome;
	mado <- mean(abs(residual));
	totresidual <- totresidual/contresidual+abs(residual);
	rnames <- rownames(data);

	if ((removeOutliers>0) && (type=="LM") && (binoutcome==FALSE ))
	{
		made <- mean(abs(totresidual),trim = 0.10); # trim @ 10% to remove outliers from estimation of mean MAD
#		cat("Bagged MAD: ", made,"\n");
		gooddatapoints <- abs(totresidual) < removeOutliers*made;
		if (sum(gooddatapoints)!=nrow(data))
		{
			reducededset <- data[gooddatapoints,];
			upmodel <- baggedModel(modelFormulas,reducededset,type,Outcome,timeOutcome,pvalue,backElimination,frequencyThreshold,removeOutliers=0);
			if (!any(is.na(upmodel$bagged.model$coefficients)))
			{
				if (length(model$coefficients) == length(upmodel$bagged.model$coefficients))
				{
					model$coefficients <- upmodel$bagged.model$coefficients;
				}
				else
				{
					for (n in 1:length(upmodel$bagged.model$coefficients))
					{
						model$coefficients[n] <- upmodel$bagged.model$coefficients[n];
					}
				}
				frma <- upmodel$formula;
				zvalues <- upmodel$zvalues;
				baggedPredict <- predictForFresa(model,data,predtype);
				residual <- baggedPredict-theoutcome;
				made <- mean(abs(residual));
				if (!any(is.na(residual)))
				{
					residual <- abs(residual/made);
					residual <- residual*(residual<10)+10*(residual>10);
					plot(residual~theoutcome,main="Outliers");
					for (i in 1:length(gooddatapoints))
					{
						if (!gooddatapoints[i])
						{
							text(theoutcome[i],residual[i],rnames[i],pos=1,cex=0.7);
						}
					}
				}
				else
				{
					cat("Warning NA residual in Bagging\n");
					print(residual);
					print(summary(model))
				}
				cat("Start MAD:", mado," Reduced MAD:",upmodel$MAD," Final MAD:",made,", Removed: ",rownames(data[!gooddatapoints,]),"\n")
				mado <- made;
			}
			else
			{
				cat("Warning NA coefficents in Bagging\n");
			}
		}
	}
	
  	result <- list(bagged.model=model,
				   formula=frma,
				   frequencyTable=VarFrequencyTable,
				   averageSize=avgsize,
				   zvalues=zvalues,
				   reducedDataSet= reducededset,
				   MAD=mado
				   );
  
	return (result);
}
