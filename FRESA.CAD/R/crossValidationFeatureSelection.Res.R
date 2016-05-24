crossValidationFeatureSelection_Res <-
function(size=10,fraction=1.0,pvalue=0.05,loops=100,covariates="1",Outcome,timeOutcome="Time",variableList,data,maxTrainModelSize=10,type=c("LM","LOGIT","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),loop.threshold=10,startOffset=0,elimination.bootstrap.steps=25,trainFraction=0.67,trainRepetition=9,elimination.pValue=0.05,setIntersect=1,interaction=c(1,1),update.pvalue=c(0.05,0.05),unirank=NULL,print=TRUE,plots=TRUE,zbaggRemoveOutliers=4.0)
{

if (!requireNamespace("cvTools", quietly = TRUE)) {
   install.packages("cvTools", dependencies = TRUE)
} 

if (!requireNamespace("glmnet", quietly = TRUE)) {
   install.packages("glmnet", dependencies = TRUE)
} 


	enetSamples <- NULL;
	enetTrainSamples <- NULL;
	totSamples <- NULL;
	totTrainSamples <- NULL;
	Full.totTrainSamples <- NULL;
	uniTrainMSS <- NULL;
	uniTestMSS <- NULL;
	
	K <- as.integer(1.0/(1.0-trainFraction) + 0.5);
	mOrderSel = interaction[1];
	mOrderUpdate = mOrderSel;
	if (length(interaction)>1) 
	{
		mOrderUpdate=interaction[2];
	}

#	cat(type,"\n")

	Fullsammples <- nrow(data);
	if ( K > Fullsammples) K=Fullsammples


	shortVarList <- as.vector(variableList[1:size,1]);
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
	shortVarList <- as.vector(rownames(table(varlist)))
	if (type=="LM")
	{
		Fullenet <- try(glmnet::cv.glmnet(as.matrix(data[,shortVarList]),as.vector(data[,Outcome]),family="gaussian"));
	}
	else
	{
		Fullenet <- try(glmnet::cv.glmnet(as.matrix(data[,shortVarList]),as.vector(data[,Outcome]),family="binomial"));
	}
	if (inherits(Fullenet, "try-error"))
	{
		cat("enet Error")
		Fullenet <- NULL;
	}
	else
	{
		cenet <- as.matrix(coef(Fullenet))
		print(LASSOVariables <- list(names(cenet[as.vector(cenet[,1]>0),])))
	}


	Full_CurModel_S <- ForwardSelection.Model.Res(size=size,fraction=fraction,pvalue=pvalue,loops=loops,covariates=covariates,Outcome=Outcome,variableList=variableList,data=data,maxTrainModelSize=maxTrainModelSize,type=type,testType=testType,timeOutcome=timeOutcome,loop.threshold=loop.threshold,interaction=mOrderSel)
	bagg <- baggedModel(Full_CurModel_S$formula.list,data,type,Outcome,timeOutcome,removeOutliers=zbaggRemoveOutliers);
	nTrainSet <- bagg$reducedDataSet;	
	if (nrow(data)!=nrow(nTrainSet)) 
	{
		Full_CurModel_S <- ForwardSelection.Model.Res(size=size,fraction=fraction,pvalue=pvalue,loops=loops,covariates=covariates,Outcome=Outcome,variableList=variableList,data=nTrainSet,maxTrainModelSize=maxTrainModelSize,type=type,testType=testType,timeOutcome=timeOutcome,loop.threshold=loop.threshold,interaction=mOrderSel)
	}


	Full_UCurModel_S <- NULL;
	Full_redCurmodel_S <- NULL;
	if (length(Full_CurModel_S$var.names)==0)
	{
		stop("no model found\n");
	}


	if (loops>1)
	{
		Full_UCurModel_S <- updateModel.Res(Outcome=Outcome,covariates=covariates,pvalue=update.pvalue,VarFrequencyTable=Full_CurModel_S$ranked.var,variableList=variableList,data=nTrainSet,type=type,testType=testType,timeOutcome=timeOutcome,interaction=mOrderUpdate,bootLoops=elimination.bootstrap.steps)
	}
	else
	{
		Full_UCurModel_S <- updateModel.Res(Outcome=Outcome,covariates=covariates,pvalue=update.pvalue,VarFrequencyTable=Full_CurModel_S$ranked.var,variableList=variableList,data=nTrainSet,type=type,testType=testType,timeOutcome=timeOutcome,interaction=mOrderUpdate,bootLoops=elimination.bootstrap.steps)
#		Full_UCurModel_S <- Full_CurModel_S;
	}



	modsize <- length(as.list(attr(terms(Full_UCurModel_S$formula),'term.labels')));
	adjsize <- min(mOrderUpdate*Full_CurModel_S$average.formula.size/pvalue,ncol(data));
	if (adjsize<2) adjsize=2;
	if (elimination.bootstrap.steps < 2 )
	{
		Full_redCurmodel_S <- backVarElimination_Res(object=Full_UCurModel_S$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=nTrainSet,startOffset=startOffset,type=type,testType=testType,setIntersect=setIntersect,adjsize=adjsize);
		Full_redCurmodel_S$bootCV  <- bootstrapValidation_Res(fraction,100,Full_redCurmodel_S$back.formula,Outcome,nTrainSet,type,plots=plots)			
	}
	else
	{		
#		cat("In Reduction \n")
		Full_redCurmodel_S <- bootstrapVarElimination_Res(object=Full_UCurModel_S$final.model,pvalue=elimination.pValue,Outcome=Outcome,
													data=nTrainSet,startOffset=startOffset,type=type,
													testType=testType,loops=elimination.bootstrap.steps,setIntersect=setIntersect,print=print,plots=plots,adjsize=adjsize);
	}

	cat ("Before BH:",as.character(Full_redCurmodel_S$beforeFSC.formula)[3],"\n")
	cat ("B:SWiMS  :",as.character(Full_redCurmodel_S$back.formula)[3],"\n")
	
	
#	cat(format(Full_redCurmodel_S$back.formula)," <-Back formula\n");
#	print(summary(Full_redCurmodel_S$back.model));
	
	
	formulas <- vector();
	BeforeBHFormulas <- vector();
	ForwardFormulas <- vector();
	baggFormulas <- vector();

	vtrainRMS <- vector();
	vblindRMS <- vector();

	vtrainSpearman <- vector();
	vtrainPearson <- vector();

	FullvtrainRMS <- vector();
	FullvblindRMS <- vector();
	FullvtrainSpearman <- vector();
	FullvtrainPearson <- vector();

	blindFoldPearson <-  vector();
	blindFoldSpearman <- vector();
	blindFoldCstat <- vector();
	blindFoldMS <- vector();
	CVBlindPearson  <- vector();
	CVBlindSpearman  <- vector();
	CVBlindRMS <- vector();

	inserted = 0;
	lastmodel = 0;

	for (i in 1:trainRepetition)
	{
		j <- 1 + ((i-1) %% K)
		if ( j == 1)
		{
#			lowFolds <- cvTools::cvFolds(nrow(lowsample), K, type = "random");
#			highFolds <- cvTools::cvFolds(nrow(hihgsample), K, type = "random");
#			lowFolds <- cvTools::cvFolds(nrow(lowsample), K,1,  "random");
#			highFolds <- cvTools::cvFolds(nrow(hihgsample), K,1, "random");
			sampleFolds <- cvTools::cvFolds(nrow(data), K,1, "random");
		}
		
#		lowTrainSet <- lowsample[lowFolds$subsets[lowFolds$which != j,],];
#		highTainSet <- hihgsample[highFolds$subsets[highFolds$which != j,],];
#		lowTestSet <- lowsample[lowFolds$subsets[lowFolds$which == j,],];
#		highTestSet <- hihgsample[highFolds$subsets[highFolds$which == j,],];
		
#		TrainSet <- rbind(lowTrainSet,highTainSet);
#		BlindSet <- rbind(lowTestSet,highTestSet);
		

		TrainSet <- data[sampleFolds$subsets[sampleFolds$which != j,],];
		BlindSet <- data[sampleFolds$subsets[sampleFolds$which == j,],];
		
		blindsampleidx <- as.vector(rownames(BlindSet));
		sampleidx <- as.vector(rownames(TrainSet));

		if (!is.null(unirank))
		{
			variableList <- update.uniRankVar(unirank,data=TrainSet,FullAnalysis=FALSE)$orderframe;
			shortVarList <- as.vector(variableList[1:size,1]);
#			print(shortVarList)
			varlist <- vector();
			for (nn in 1:length(shortVarList))
			{
				varlist <- append(varlist,str_replace_all(unlist(strsplit(
								str_replace_all(
									str_replace_all(
										str_replace_all(
											str_replace_all(shortVarList[nn],"I\\("," ")
										,"\\("," ")
									,">","\\*")
								,"<","\\*")
						,"\\*"))[1]," ",""))
			}
			shortVarList <- as.vector(rownames(table(varlist)))
#			print(shortVarList)
		}




		cat("Samples Train :",nrow(TrainSet),"Samples Test :",nrow(BlindSet),"\n");
		cat ("Loop :",i,"\n")

		if (!is.null(Fullenet))
		{
			if (type=="LM")
			{
				foldenet <- try(glmnet::cv.glmnet(as.matrix(TrainSet[,shortVarList]),as.vector(TrainSet[,Outcome]),family="gaussian"));
			}
			else
			{
				foldenet <- try(glmnet::cv.glmnet(as.matrix(TrainSet[,shortVarList]),as.vector(TrainSet[,Outcome]),family="binomial"));
			}
			cenet <- as.matrix(coef(foldenet))
			LASSOVariables[[i+1]] <- names(cenet[as.vector(cenet[,1]>0),])
			if (i == 1)
			{
				enetSamples <- cbind(BlindSet[,Outcome],predict(foldenet,as.matrix(BlindSet[,shortVarList])),i);
				enetTrainSamples <- cbind(TrainSet[,Outcome],predict(foldenet,as.matrix(TrainSet[,shortVarList])),i);
			}
			else
			{
				enetSamples <- rbind(enetSamples,cbind(BlindSet[,Outcome],predict(foldenet,as.matrix(BlindSet[,shortVarList])),i));
				enetTrainSamples <- rbind(enetTrainSamples,cbind(TrainSet[,Outcome],predict(foldenet,as.matrix(TrainSet[,shortVarList])),i));
			}
##			print(LASSOVariables)
		}

		
		
		par(mfrow=c(1,1))

#		cat(type,"\n")

#		cat(Full_redCurmodel_S$back.formula," <-Back formula\n");

		CurModel_S <- ForwardSelection.Model.Res(size=size,fraction=fraction,pvalue=pvalue,loops=loops,covariates=covariates,Outcome=Outcome,variableList=variableList,data=TrainSet,maxTrainModelSize=maxTrainModelSize,type=type,testType=testType,timeOutcome=timeOutcome,loop.threshold=loop.threshold,interaction=mOrderSel)
		if (length(CurModel_S$var.names)>0)
		{
			bagg <- baggedModel(CurModel_S$formula.list,TrainSet,type,Outcome,timeOutcome,removeOutliers=zbaggRemoveOutliers);
			nTrainSet <- bagg$reducedDataSet;
#			cat("Train Size: ",nrow(TrainSet)," Reduced Train Size:",nrow(nTrainSet),"\n")
			if (nrow(TrainSet)!=nrow(nTrainSet)) 
			{
				cat("Before Outilier Initial Forward Model:",as.character(CurModel_S$formula)[3],"\n");
				CurModel_S <- ForwardSelection.Model.Res(size=size,fraction=fraction,pvalue=pvalue,loops=loops,covariates=covariates,Outcome=Outcome,variableList=variableList,data=nTrainSet,maxTrainModelSize=maxTrainModelSize,type=type,testType=testType,timeOutcome=timeOutcome,loop.threshold=loop.threshold,interaction=mOrderSel)
				bagg <- baggedModel(CurModel_S$formula.list,nTrainSet,type,Outcome,timeOutcome,removeOutliers=0);
				cat("After Outlier Removing Forward Model :",as.character(CurModel_S$formula)[3],"\n");
			}
			inserted = inserted +1;
			if (loops>1)
			{
				UCurModel_S <- updateModel.Res(Outcome=Outcome,covariates=covariates,pvalue=update.pvalue,VarFrequencyTable=CurModel_S$ranked.var,variableList=variableList,data=nTrainSet,type=type,testType=testType,timeOutcome=timeOutcome,interaction=mOrderUpdate,bootLoops=elimination.bootstrap.steps)
			}
			else
			{
#				UCurModel_S <- CurModel_S;
				UCurModel_S <- updateModel.Res(Outcome=Outcome,covariates=covariates,pvalue=update.pvalue,VarFrequencyTable=CurModel_S$ranked.var,variableList=variableList,data=nTrainSet,type=type,testType=testType,timeOutcome=timeOutcome,interaction=mOrderUpdate,bootLoops=elimination.bootstrap.steps)
			}
			modsize <- length(as.list(attr(terms(UCurModel_S$formula),'term.labels')));
			adjsize <- min(mOrderUpdate*CurModel_S$average.formula.size/pvalue,ncol(data));
			if (adjsize<2) adjsize=2;
			if (elimination.bootstrap.steps < 2 )
			{
				redCurmodel_S <- backVarElimination_Res(object=UCurModel_S$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=nTrainSet,startOffset=startOffset,type=type,testType=testType,setIntersect=setIntersect,adjsize=adjsize);
				redCurmodel_S$bootCV  <- bootstrapValidation_Res(fraction,100,redCurmodel_S$back.formula,Outcome,nTrainSet,type,plots=plots)
			}
			else
			{
				redCurmodel_S <- bootstrapVarElimination_Res(object=UCurModel_S$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=nTrainSet,startOffset=startOffset,type=type,testType=testType,loops=elimination.bootstrap.steps,setIntersect=setIntersect,print=print,plots=plots,adjsize=adjsize);
			}

			Full_model <- modelFitting(Full_redCurmodel_S$back.model,TrainSet,type)


			redfoldmodel.BBH <- redCurmodel_S$beforeFSC.model;
			redfoldmodel <- redCurmodel_S$bootCV$boot.model;
#			redfoldmodel.BBH <- redCurmodel_S$beforeFSC.model;
#			redfoldmodel <- redCurmodel_S$back.model;
			
			cat ("Update   :",as.character(UCurModel_S$formula)[3],"\n")
			cat ("Before BH:",as.character(redCurmodel_S$beforeFSC.formula)[3],"\n")
			cat ("B:SWiMS  :",as.character(redCurmodel_S$back.formula)[3],"\n")
			
			if (!is.null(redfoldmodel))
			{
			
				predictTest <- predictForFresa(redfoldmodel,BlindSet,"linear");
				predictTest.BBH <- predictForFresa(redfoldmodel.BBH,BlindSet,"linear");
				predictTest.ForwardModel <- predictForFresa(UCurModel_S$final.model,BlindSet,"linear");
				Full_predictTest <- predictForFresa(Full_model,BlindSet,"linear");
				predictTrain <- predictForFresa(redfoldmodel,TrainSet,"linear");
				Full_predictTrain <- predictForFresa(Full_model,TrainSet,"linear");
				
	#			bagg <- baggedModel(CurModel_S$formula.list,TrainSet,"LM",Outcome); 
				baggedForwardPredict <- predictForFresa(bagg$bagged.model,BlindSet,"linear");
				medianPred <- medianPredict(CurModel_S$formula.list,nTrainSet,BlindSet, predictType = "linear",type = type)$medianPredict



				trainResiduals <- residualForFRESA(redfoldmodel,TrainSet,Outcome);
				blindResiduals <- residualForFRESA(redfoldmodel,BlindSet,Outcome);
				FulltrainResiduals <- residualForFRESA(Full_model,TrainSet,Outcome);
				FullblindResiduals <- residualForFRESA(Full_model,BlindSet,Outcome);


				if (inserted == 1)
				{
					totSamples <- cbind(BlindSet[,Outcome],predictTest,i,blindResiduals,medianPred,baggedForwardPredict,predictTest.ForwardModel,predictTest.BBH);
					rownames(totSamples) <- blindsampleidx;
					Full.totSamples <- cbind(BlindSet[,Outcome],Full_predictTest,i,FullblindResiduals);
					rownames(Full.totSamples) <- blindsampleidx;
					totTrainSamples <- cbind(TrainSet[,Outcome],predictTrain,i,trainResiduals);
					rownames(totTrainSamples) <- sampleidx;
					Full.totTrainSamples <- cbind(TrainSet[,Outcome],Full_predictTrain,i,FulltrainResiduals);
					rownames(Full.totTrainSamples) <- sampleidx;
				}
				else
				{
					px <- cbind(BlindSet[,Outcome],predictTest,i,blindResiduals,medianPred,baggedForwardPredict,predictTest.ForwardModel,predictTest.BBH);
					rownames(px) <- blindsampleidx;
					totSamples <- rbind(totSamples,px);
					px <- cbind(BlindSet[,Outcome],Full_predictTest,i,FullblindResiduals);
					rownames(px) <- blindsampleidx;
					Full.totSamples <- rbind(Full.totSamples,px);
					px <- cbind(TrainSet[,Outcome],predictTrain,i,trainResiduals);
					rownames(px) <- sampleidx;
					totTrainSamples <- rbind(totTrainSamples,px);
					px <- cbind(TrainSet[,Outcome],Full_predictTrain,i,FulltrainResiduals);
					rownames(px) <- sampleidx;
					Full.totTrainSamples <- rbind(Full.totTrainSamples,px);

				}
				formulas <- append(formulas,as.character(redCurmodel_S$back.formula)[3]);
				BeforeBHFormulas <- append(BeforeBHFormulas,as.character(redCurmodel_S$beforeFSC.formula)[3]);
				ForwardFormulas <- append(ForwardFormulas,as.character(UCurModel_S$formula)[3]);
				baggFormulas <- append(baggFormulas,bagg$formula);
					
				
				trainRMS <- sqrt(sum(trainResiduals^2)/nrow(TrainSet));
				trainPearson <- cor.test(TrainSet[,Outcome], predictTrain, method = "pearson",na.action=na.omit,exact=FALSE)$estimate
				trainSpearman <- cor.test(TrainSet[,Outcome], predictTrain, method = "spearman",na.action=na.omit,exact=FALSE)$estimate

				FulltrainRMS <- sqrt(sum(FulltrainResiduals^2)/nrow(TrainSet));
				FulltrainPearson <- cor.test(TrainSet[,Outcome], Full_predictTrain, method = "pearson",na.action=na.omit,exact=FALSE)$estimate
				FulltrainSpearman <- cor.test(TrainSet[,Outcome], Full_predictTrain, method = "spearman",na.action=na.omit,exact=FALSE)$estimate

				blindRMS <- sqrt(sum(blindResiduals^2)/nrow(BlindSet));
				FullblindRMS <- sqrt(sum(FullblindResiduals^2)/nrow(BlindSet));

				if (nrow(BlindSet)>5)
				{
					foldPearson <- cor.test(BlindSet[,Outcome], predictTest, method = "pearson",na.action=na.omit,exact=FALSE)$estimate
					foldSpearman <- cor.test(BlindSet[,Outcome], predictTest, method = "spearman",na.action=na.omit,exact=FALSE)$estimate				
					cstat <- rcorr.cens(predictTest,BlindSet[,Outcome], outx=FALSE)[1];
					foldRMS <- sum(blindResiduals^2)/(nrow(BlindSet)-1);
					cat("Fold RMS: ",sqrt(foldRMS),"Fold Test Pearson: ", foldPearson, "Fold Test Spearman: ",foldSpearman,"Fold Cstat:",cstat,"\n");

					blindFoldMS <- append(blindFoldMS,foldRMS);
					blindFoldPearson <- append(blindFoldPearson,foldPearson);
					blindFoldSpearman <- append(blindFoldSpearman,foldSpearman);
					blindFoldCstat <- append(blindFoldCstat,cstat);
	#				cat("Accu RMS: ",sqrt(mean(blindFoldMS)),"Accu Test Pearson: ", mean(blindFoldPearson), "Accu Test Spearman: ",mean(blindFoldSpearman),"Accu Cstat:",mean(blindFoldCstat),"\n");
				}
	#			cat("After RMS\n");

	# univariate analysis of top model residuals
				uniEval <- getVar.Res(Full_redCurmodel_S$back.formula,TrainSet,Outcome,type = type,testData=BlindSet);
				if (i==1)
				{
					uniTrainMSS <- rbind(uniEval$unitrainMSS);
					uniTestMSS <- rbind(uniEval$unitestMSS);
				}
				else
				{
					uniTrainMSS <- rbind(uniTrainMSS,uniEval$unitrainMSS);
					uniTestMSS <- rbind(uniTestMSS,uniEval$unitestMSS);
				}


				if (nrow(totSamples)>5)
				{
					blindPearson <- cor.test(totSamples[,1], totSamples[,2], method = "pearson",na.action=na.omit,exact=FALSE)$estimate
					blindSpearman <- cor.test(totSamples[,1], totSamples[,2], method = "spearman",na.action=na.omit,exact=FALSE)
					blindForwardSpearman <- cor.test(totSamples[,1], totSamples[,6], method = "spearman",na.action=na.omit,exact=FALSE)$estimate

					cstat <- rcorr.cens(totSamples[,2],totSamples[,1], outx=FALSE)[1];

					FullblindPearson <- cor.test(Full.totSamples[,1], Full.totSamples[,2], method = "pearson",na.action=na.omit,exact=FALSE)$estimate
					FullblindSpearman <- cor.test(Full.totSamples[,1], Full.totSamples[,2], method = "spearman",na.action=na.omit,exact=FALSE)

					AcumRMS <- sqrt(sum(totSamples[,4]^2)/nrow(totSamples));
					AcumFullRMS <- sqrt(sum(Full.totSamples[,4]^2)/nrow(totSamples));


					cat("Samples: ", nrow(totSamples),"Full RMS:",AcumFullRMS," Accumulated Blind RMS: ", AcumRMS," c-index : ",cstat,"\n");
					cat("Full Blind RMS: ", FullblindRMS, " Full Train RMS: ",FulltrainRMS,"\n");
					cat("Blind Pearson: ", blindPearson, " Train Pearson: ",trainPearson,"\n");
					cat("Blind Spearman: ", blindSpearman$estimate,"(", blindSpearman$p.value,")  Train Spearman: ",trainSpearman,"Forward Bagged Spearman:",blindForwardSpearman,"\n");
					cat("Full Blind Pearson: ", FullblindPearson , " Full Train Pearson: ",FulltrainPearson,"\n");
					cat("Full Blind Spearman: ", FullblindSpearman$estimate, "(",FullblindSpearman$p.value,") Full Train Spearman: ",FulltrainSpearman,"\n");

					

				}
				
				vblindRMS <- append(vblindRMS,blindRMS);
				FullvblindRMS <- append(FullvblindRMS,FullblindRMS);
				
				vtrainRMS <- append(vtrainRMS,trainRMS);
				vtrainSpearman <- append(vtrainSpearman,trainSpearman);
				vtrainPearson <- append(vtrainPearson,trainPearson);

				FullvtrainRMS <- append(FullvtrainRMS,FulltrainRMS);
				FullvtrainSpearman <- append(FullvtrainSpearman,FulltrainSpearman);
				FullvtrainPearson <- append(FullvtrainPearson,FulltrainPearson);
			}


		}
		if ( (i %% K) == 0)
		{
			foldtest <- totSamples[totSamples[,3]>lastmodel,];
			CVBlindPearson <- append(CVBlindPearson,cor.test(foldtest[,1], foldtest[,2], method = "pearson",na.action=na.omit,exact=FALSE)$estimate);
			CVBlindSpearman <- append(CVBlindSpearman,cor.test(foldtest[,1], foldtest[,2], method = "spearman",na.action=na.omit,exact=FALSE)$estimate);
			CVBlindRMS <- append(CVBlindRMS,sqrt(sum((foldtest[,1]-foldtest[,2])^2)/nrow(foldtest)));
			lastmodel = i;
		}
	}
			
#	print(LASSOVariables)

	colnames(totSamples) <- c("Outcome","Prediction","Model","Residuals","Median","Bagged","Forward","Backwards");
	totSamples <- as.data.frame(totSamples);

	colnames(Full.totSamples) <- c("Outcome","Prediction","Model","Residuals");
	Full.totSamples <- as.data.frame(Full.totSamples);

	colnames(totTrainSamples) <- c("Outcome","Prediction","Model","Residuals");
	totTrainSamples <- as.data.frame(totTrainSamples);

	colnames(Full.totTrainSamples) <- c("Outcome","Prediction","Model","Residuals");
	Full.totTrainSamples <- as.data.frame(Full.totTrainSamples);

	
	BSWiMS.ensemble.prediction <- NULL
	bsta <- boxplot(totSamples$Prediction~rownames(totSamples),plot=FALSE)
	sta <- cbind(bsta$stats[3,])
	rownames(sta) <- bsta$names
#	BSWiMS.ensemble.prediction <- cbind(data[,Outcome],sta[rownames(data),])
	BSWiMS.ensemble.prediction <- cbind(data[rownames(sta),Outcome],sta)
	colnames(BSWiMS.ensemble.prediction) <- c("Outcome","Prediction");
	BSWiMS.ensemble.prediction <- as.data.frame(BSWiMS.ensemble.prediction);


	if (!is.null(enetSamples))
	{
		colnames(enetSamples) <- c("Outcome","Prediction","Model");
		enetSamples <- as.data.frame(enetSamples);
		colnames(enetTrainSamples) <- c("Outcome","Prediction","Model");
		enetTrainSamples <- as.data.frame(enetTrainSamples);
	}

	if (!is.null(uniTrainMSS))
	{
		uniTrainMSS <- as.data.frame(uniTrainMSS);
		uniTestMSS <- as.data.frame(uniTestMSS);
		colnames(uniTrainMSS) <-  attr(terms(Full_redCurmodel_S$back.formula),'term.labels');
		colnames(uniTestMSS) <-  attr(terms(Full_redCurmodel_S$back.formula),'term.labels');
	}

	
	blindRMS <- sqrt(sum((totSamples$Residuals)^2)/nrow(totSamples));
	blindPearson <- cor.test(totSamples$Outcome, totSamples$Prediction, method = "pearson",na.action=na.omit,exact=FALSE)$estimate
	blindSpearman <- cor.test(totSamples$Outcome, totSamples$Prediction, method = "spearman",na.action=na.omit,exact=FALSE)$estimate

	FullblindRMS <- sqrt(sum((Full.totSamples$Residuals)^2)/nrow(totSamples));
	FullblindPearson <- cor.test(Full.totSamples$Outcome, Full.totSamples$Prediction, method = "pearson",na.action=na.omit,exact=FALSE)$estimate
	FullblindSpearman <- cor.test(Full.totSamples$Outcome, Full.totSamples$Prediction, method = "spearman",na.action=na.omit,exact=FALSE)$estimate

	cstat <- rcorr.cens(totSamples$Prediction,totSamples$Outcome, outx=FALSE)[1];

	cat("##### CV of Initial Model (Biased) ###### \n");
	cat("Full Blind RMS: ", FullblindRMS,"\n")
	cat("Full Blind Spearman: ", FullblindSpearman,"\n")
	cat("Full Blind Pearson: ", FullblindPearson,"\n")
	if (length(blindFoldMS)>0)
	{
		cat("##### By Fold Analysis ###### \n Samples: ",length(blindFoldMS),"\n");
		cat("Sampled RMS: ", sqrt(mean(blindFoldMS)),"\n")
		cat("Mean Blind Spearman: ", mean(blindFoldSpearman)," (",sd(blindFoldSpearman),")\n")
		cat("Mean Blind Pearson: ", mean(blindFoldPearson)," (",sd(blindFoldPearson),")\n")
		cat("Mean Blind cstat: ",mean(blindFoldCstat)," (",sd(blindFoldCstat),")\n")
	}
	cat("##### Full Set Coherence Analysis ###### \n");
	cat("Blind  RMS: ", blindRMS," c-index : ",cstat,"\n");
	cat("Blind  Spearman: ", blindSpearman,"\n")
	cat("Blind  Pearson: ", blindPearson,"\n")
	if (length(CVBlindPearson)>0)
	{
		cat("##### Models Coherence Analysis ###### \n Samples: ",length(CVBlindPearson),"\n");
		cat("Mean Blind Spearman: ", mean(CVBlindSpearman)," (",sd(CVBlindSpearman),")\n")
		cat("Mean Blind Pearson: ", mean(CVBlindPearson)," (",sd(CVBlindPearson),")\n")
	}

	
	
	result <- list(		formula.list=formulas, 
						Models.testPrediction=totSamples,
						FullBSWiMS.testPrediction=Full.totSamples,
						BSWiMS=Full_redCurmodel_S,
						forwardSelection=Full_CurModel_S,
						updatedforwardModel=Full_UCurModel_S,
						testRMSE = blindRMS,
						testPearson = blindPearson,
						testSpearman = blindSpearman,
						FullTestRMSE = FullblindRMS,
						FullTestPearson = FullblindPearson,
						FullTestSpearman = FullblindSpearman,
						trainRMSE = vtrainRMS,
						trainPearson = vtrainPearson,
						trainSpearman = vtrainSpearman,
						FullTrainRMS = FullvtrainRMS,
						FullTrainPearson = FullvtrainPearson,
						FullTrainSpearman = FullvtrainSpearman,
						byFoldTestMS = blindFoldMS,
						byFoldTestSpearman = blindFoldSpearman,
						byFoldTestPearson = blindFoldPearson,
						byFoldCstat = blindFoldCstat,
						testRMSEAtFold = vblindRMS,
						FullTestRMSEAtFold = FullvblindRMS,
						Fullenet=Fullenet,
						LASSO.testPredictions=enetSamples,
						LASSOVariables=LASSOVariables,
						CVBlindPearson=CVBlindPearson,
						CVBlindSpearman=CVBlindSpearman,
						CVBlindRMS=CVBlindRMS,
						Models.trainPrediction=totTrainSamples,
						FullBSWiMS.trainPrediction=Full.totTrainSamples,
						LASSO.trainPredictions=enetTrainSamples,
						uniTrainMSS = uniTrainMSS,
						uniTestMSS = uniTestMSS,
						BSWiMS.ensemble.prediction = BSWiMS.ensemble.prediction,
						BeforeBHFormulas.list = BeforeBHFormulas,
						ForwardFormulas.list = ForwardFormulas,
						baggFormulas.list = baggFormulas

					);
	return (result)
}
