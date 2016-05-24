crossValidationFeatureSelection_Bin <-
function(size=10,fraction=1.0,pvalue=0.05,loops=100,covariates="1",Outcome,timeOutcome="Time",variableList,data,maxTrainModelSize=10,type=c("LM","LOGIT","COX"),selectionType=c("zIDI","zNRI","Both"),loop.threshold=10,startOffset=0,elimination.bootstrap.steps=25,trainFraction=0.67,trainRepetition=9,elimination.pValue=0.05,CVfolds=10,bootstrap.steps=25,interaction=c(1,1),nk=0,unirank=NULL,print=TRUE,plots=TRUE)
{

if (!requireNamespace("cvTools", quietly = TRUE)) {
   install.packages("cvTools", dependencies = TRUE)
} 

if (!requireNamespace("glmnet", quietly = TRUE)) {
   install.packages("glmnet", dependencies = TRUE)
} 


	enetSamples <- NULL;
	enetTrainSamples <- NULL;

	casesample = subset(data,get(Outcome)  == 1);
	controlsample = subset(data,get(Outcome) == 0);

	casesamplesize <- nrow(casesample);

	controlsamplesize <- nrow(controlsample);
	

	K <- as.integer(1.0/(1.0-trainFraction) + 0.5);


	acc = 0.0;
	sen = 0.0;
	spe = 0.0;
	sizecases = 0;
	sizecontrol = 0;
	totsize = 0;
	paracc = 0;
	psen = 0;
	pspe = 0;

	Full.acc = 0.0;
	Full.sen = 0.0;
	Full.spe = 0.0;
	Full.paracc = 0;
	Full.psen = 0;
	Full.pspe = 0;



	formulas <- vector();
	BeforeBHFormulas <- vector();
	ForwardFormulas <- vector();
	baggFormulas <- vector();
	trainCorrelations <- vector();
	trainAccuracy <- vector();
	trainSensitivity <- vector();
	trainSpecificity <- vector();
	trainAUC <- vector();
	testAccuracy <- vector();
	testSensitivity <- vector();
	testSpecificity <- vector();
	testAUC <- vector();

	blindCorrelations <- vector();
	WholeFoldBlindAccuracy <- vector();
	WholeFoldBlindSpecificity <- vector();
	WholeFoldBlindSensitivity <- vector();
	WholeFoldBlindAUC <- vector();
	FoldBlindAccuracy <- vector();
	FoldBlindSpecificity <- vector();
	FoldBlindSensitivity <- vector();
	TopUniCoherenceTest <- vector();
	selection.pValue <- pvalue;
#	update.pValue <- c(pvalue*pvalue,pvalue);
	update.pValue <- c(pvalue,0.99*pvalue);

	CVselection.pValue <- pvalue;
	CVelimination.pValue <- elimination.pValue;
	CVupdate.pValue <- update.pValue


	par(mfrow=c(1,1))
	
	mOrderSel = interaction[1];
	mOrderUpdate = mOrderSel;
	if (length(interaction)>1) mOrderUpdate=interaction[2]
	
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
#	print(shortVarList)
	Fullenet <- try(glmnet::cv.glmnet(as.matrix(data[,shortVarList]),as.vector(data[,Outcome]),family="binomial"));
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
	

	
	
	selType = selectionType;
	if (selectionType=="Both") selType ="zIDI";

	CurModel_Full <- ForwardSelection.Model.Bin(size=size,fraction=fraction,pvalue=CVselection.pValue,loops=loops,covariates=covariates,Outcome=Outcome,variableList=variableList,data=data,maxTrainModelSize=maxTrainModelSize,type=type,timeOutcome=timeOutcome,selectionType=selectionType,loop.threshold=loop.threshold,interaction=mOrderSel)
	
	if (length(CurModel_Full$var.names)==0)
	{
		stop("no initial model found\n");
	}

	if (loops>1)
	{
		UCurModel_Full <- updateModel.Bin(Outcome=Outcome,covariates=covariates,pvalue=CVupdate.pValue,VarFrequencyTable=CurModel_Full$ranked.var,variableList=variableList,data=data,type=type,lastTopVariable= 0,timeOutcome=timeOutcome,selectionType=selType,interaction=mOrderUpdate,numberOfModels=0,bootLoops=elimination.bootstrap.steps)
	}
	else
	{
		UCurModel_Full <- updateModel.Bin(Outcome=Outcome,covariates=covariates,pvalue=CVupdate.pValue,VarFrequencyTable=CurModel_Full$ranked.var,variableList=variableList,data=data,type=type,lastTopVariable= 0,timeOutcome=timeOutcome,selectionType=selType,interaction=mOrderUpdate,numberOfModels=0,bootLoops=elimination.bootstrap.steps)
#		UCurModel_Full <- CurModel_Full;
	}
	# the estimated number of independent features
	modsize <- length(as.list(attr(terms(UCurModel_Full$formula),'term.labels')));
	adjsize <- min(mOrderUpdate*CurModel_Full$average.formula.size/pvalue,ncol(data));
	if (adjsize<2) adjsize=2;
	if (elimination.bootstrap.steps>1)
	{
		redCurmodel_Full <- bootstrapVarElimination_Bin(object=UCurModel_Full$final.model,pvalue=CVelimination.pValue,Outcome=Outcome,data=data,startOffset=startOffset,type=type,selectionType=selType,loops=elimination.bootstrap.steps,fraction=fraction,print=print,plots=plots,adjsize=adjsize);
	}
	else
	{
		
		redCurmodel_Full <- backVarElimination_Bin(object=UCurModel_Full$final.model,pvalue=CVelimination.pValue,Outcome=Outcome,data=data,startOffset=startOffset,type=type,selectionType=selType,adjsize=adjsize);
	}
	
	Full_formula <- redCurmodel_Full$back.model;
	FullBootCross <- bootstrapValidation_Bin(1.0000,bootstrap.steps,Full_formula,Outcome,data,type,plots=plots)
	redBootCross <- FullBootCross;

	cat ("Update   :",as.character(UCurModel_Full$formula)[3],"\n")
	cat ("Before BH:",as.character(redCurmodel_Full$beforeFSC.formula)[3],"\n")
	cat ("B:SWiMS  :",as.character(redCurmodel_Full$back.formula)[3],"\n")

	if (is.null(FullBootCross))
	{
		stop("no initial model found\n");
	}
	if (print) summary(FullBootCross,2)

	
	inserted = 0;
	rocadded = 0;
	split.blindSen <- NULL;
	blindreboot <- NULL;
	KNNSamples <- NULL;
	Full.KNNSamples <- NULL;
	totSamples <- NULL;
	Full.totSamples <- NULL;
	totTrainSamples <- NULL;
	Full.totTrainSamples <- NULL;
	uniTrainAccuracy <- NULL;
	uniTestAccuracy <- NULL;

	Fullsammples <- min(casesamplesize,controlsamplesize);
	if ( K > Fullsammples) K=Fullsammples
	cat("Number of folds: ",K,"\n");
	specificities <- c(0.975,0.95,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.05);


	for (i in 1:trainRepetition)
	{
		j <- 1 + ((i-1) %% K)
		if ( j == 1)
		{
#			casefolds <- cvTools::cvFolds(casesamplesize, K, type = "random");
#			controlfolds <- cvTools::cvFolds(controlsamplesize, K, type = "random");
			casefolds <- cvTools::cvFolds(casesamplesize, K,1,  "random");
			controlfolds <- cvTools::cvFolds(controlsamplesize, K,1,  "random");
			cycleinsert=0;
			totalUniCor=0;
		}

		CaseTrainSet <- casesample[casefolds$subsets[casefolds$which != j,],];
		CaseBlindSet <- casesample[casefolds$subsets[casefolds$which == j,],];
		ControlTrainSet <- controlsample[controlfolds$subsets[controlfolds$which != j,],];
		ControlBlindSet <- controlsample[controlfolds$subsets[controlfolds$which == j,],];

		TrainSet <- rbind(CaseTrainSet,ControlTrainSet);
		BlindSet <- rbind(CaseBlindSet,ControlBlindSet);
		framesize <- nrow(BlindSet);
		minTrainSamples <- min(nrow(CaseTrainSet),nrow(ControlTrainSet));

		if (nk==0)
		{
			nk = 2*as.integer(sqrt(minTrainSamples/2)) + 1;
		}
			

		KnnTrainSet <- rbind(CaseTrainSet[sample(1:nrow(CaseTrainSet),minTrainSamples,replace=FALSE),],ControlTrainSet[sample(1:nrow(ControlTrainSet),minTrainSamples,replace=FALSE),])
		
		par(mfrow=c(1,1))

#		redBootCross <- bootstrapValidation_Bin(1.0000,bootstrap.steps,Full_formula,Outcome,TrainSet,type,plots=plots)
		redBootCross$boot.model <- modelFitting(Full_formula,TrainSet,type);

#		par(mfrow=c(1,1))

#		print(summary(redBootCross$boot.model))
		Full.p <- predictForFresa(redBootCross$boot.model,BlindSet, 'linear');
		
		Fullknnclass <- getKNNpredictionFromFormula(UCurModel_Full$formula,KnnTrainSet,BlindSet,Outcome,nk)


		if (!is.null(unirank))
		{
#			cat("Ranking Again \n");
			variableList <- update.uniRankVar(unirank,data=TrainSet,FullAnalysis=FALSE)$orderframe;
#			cat("Ranked Again. Size: ",nrow(variableList),":",ncol(variableList),", \n");
#			print(variableList[1:10,]);
			shortVarList <- as.vector(variableList[1:size,1]);
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
#			cat("End Ranking Again \n");


		}
		
		
		
		if (!is.null(Fullenet))
		{
#			cat("In elastic Net\n")
			foldenet <- try(glmnet::cv.glmnet(as.matrix(TrainSet[,shortVarList]),as.vector(TrainSet[,Outcome]),family="binomial"));
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
#			print(LASSOVariables)
		}




		cat ("Loop :",i,"Input Cases =",sum(data[,Outcome] > 0 ),"Input Control =",sum(data[,Outcome] == 0),"\n")
		cat ("Loop :",i,"Train Cases =",sum(TrainSet[,Outcome] > 0 ),"Train Control =",sum(TrainSet[,Outcome] == 0),"\n")
		cat ("Loop :",i,"Blind Cases =",sum(BlindSet[,Outcome] > 0 ),"Blind Control =",sum(BlindSet[,Outcome] == 0),"\n")
		cat ("K   :",nk,"KNN T Cases =",sum(KnnTrainSet[,Outcome] > 0 ),"KNN T Control =",sum(KnnTrainSet[,Outcome] == 0),"\n")

		lastinserted = inserted;
		CurModel_S <- ForwardSelection.Model.Bin(size=size,fraction=fraction,pvalue=selection.pValue,loops=loops,covariates=covariates,Outcome=Outcome,variableList=variableList,data=TrainSet,maxTrainModelSize=maxTrainModelSize,type=type,timeOutcome=timeOutcome,selectionType=selectionType,loop.threshold=loop.threshold,interaction=mOrderSel)
		if (length(CurModel_S$var.names)>0)
		{
			if (loops>1)
			{
				UCurModel_S <- updateModel.Bin(Outcome=Outcome,covariates=covariates,pvalue=update.pValue,VarFrequencyTable=CurModel_S$ranked.var,variableList=variableList,data=TrainSet,type=type,lastTopVariable= 0,timeOutcome=timeOutcome,selectionType=selType,interaction=mOrderUpdate,numberOfModels=0,bootLoops=elimination.bootstrap.steps)
			}
			else
			{
				UCurModel_S <- updateModel.Bin(Outcome=Outcome,covariates=covariates,pvalue=update.pValue,VarFrequencyTable=CurModel_S$ranked.var,variableList=variableList,data=TrainSet,type=type,lastTopVariable= 0,timeOutcome=timeOutcome,selectionType=selType,interaction=mOrderUpdate,numberOfModels=0,bootLoops=elimination.bootstrap.steps)
#				UCurModel_S <- CurModel_S;
			}
			modsize <- length(as.list(attr(terms(UCurModel_S$formula),'term.labels')));
			adjsize <- min(mOrderUpdate*CurModel_S$average.formula.size/pvalue,ncol(data));
			if (adjsize<2) adjsize=2;
			if (elimination.bootstrap.steps > 1)
			{
				redCurmodel_S <- bootstrapVarElimination_Bin(object=UCurModel_S$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=TrainSet,startOffset=startOffset,type=type,selectionType=selType,loops=elimination.bootstrap.steps,fraction=fraction,print=print,plots=plots,adjsize=adjsize);	
				redBootCross_S <- redCurmodel_S$bootCV;
			}
			else
			{
				redCurmodel_S <- backVarElimination_Bin(object=UCurModel_S$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=TrainSet,startOffset=startOffset,type=type,selectionType=selType,adjsize=adjsize);
				redBootCross_S <- bootstrapValidation_Bin(1.0000,bootstrap.steps,redCurmodel_S$back.formula,Outcome,TrainSet,type,plots=plots)
				par(mfrow=c(1,1))
			}
			 
#			redBootCross_S$boot.model <- modelFitting(redCurmodel_S$back.formula,TrainSet,type);

# lets use beforeFSC model for prediction 
			redfoldmodel.BBH <- redCurmodel_S$beforeFSC.model;
			forwardmodel <- UCurModel_S$final.model;
			redfoldmodel <- redBootCross_S$boot.model;



			cat ("Update   :",as.character(UCurModel_S$formula)[3],"\n")
			cat ("Before BH:",as.character(redCurmodel_S$beforeFSC.formula)[3],"\n")
			cat ("B:SWiMS  :",as.character(redCurmodel_S$back.formula)[3],"\n")

			if ((redCurmodel_S$lastRemoved >= 0) && !is.null(redfoldmodel))
			{

				if (print) 
				{
					cat ("\n The last CV bootstrapped model")
					s <- summary(redBootCross_S,2)
				}
			
				bagg <- baggedModel(CurModel_S$formula.list,TrainSet,type,Outcome,timeOutcome,removeOutliers=0.0); 
#				bagg <- baggedModel(CurModel_S$formula.list,TrainSet,"LM",Outcome); 
#				bagg <- baggedModel(CurModel_S$formula.list,KnnTrainSet,"LM",Outcome); 
				knnclass <- getKNNpredictionFromFormula(UCurModel_S$formula,KnnTrainSet,BlindSet,Outcome,nk)

			
				p <- predictForFresa(redfoldmodel,BlindSet, 'linear');
				p.BBH <- predictForFresa(redfoldmodel.BBH,BlindSet, 'linear');
				p.forward <- predictForFresa(forwardmodel,BlindSet, 'linear');
				baggedForwardPredict <- predictForFresa(bagg$bagged.model,BlindSet,'linear');
				medianPred <- medianPredict(CurModel_S$formula.list,TrainSet,BlindSet, predictType = "linear",type = type)$medianPredict
				
				inserted = inserted + 1
				cycleinsert = cycleinsert + 1

				tcor <- cor.test(predictForFresa(redfoldmodel,TrainSet, 'linear'),predictForFresa(redBootCross$boot.model,TrainSet, 'linear'), method = "spearman",na.action=na.omit,exact=FALSE)$estimate
				trainCorrelations <- append(trainCorrelations,tcor);
				trainAccuracy <- append(trainAccuracy,redBootCross_S$base.Accuracy);
				trainSensitivity <- append(trainSensitivity,redBootCross_S$base.Sensitivity);
				trainSpecificity <- append(trainSpecificity,redBootCross_S$base.Specificity);
				trainAUC <- append(trainAUC,mean(redBootCross_S$train.ROCAUC));


				bcor <- 0;
				if (framesize>5)
				{
					bcor <- cor.test(p, Full.p, method = "spearman",na.action=na.omit,exact=FALSE)$estimate;
					blindCorrelations <- append(blindCorrelations,bcor);
				}

				if (((sumca <- sum(BlindSet[,Outcome]>0)) > 1) && ((sumco <- sum(BlindSet[,Outcome]==0)) > 1))
				{

					atRoc <- pROC::roc(BlindSet[,Outcome], p,plot=FALSE,ci=TRUE,auc=TRUE,of='se',specificities=specificities,boot.n=100,smooth=FALSE,progress= 'none')
					splitRoc <- atRoc$ci[,2];
					FullRocBlindAUC <- pROC::roc(BlindSet[,Outcome], Full.p,plot=FALSE,auc=TRUE,ci=FALSE)$auc
					WholeFoldBlindAUC <- append(WholeFoldBlindAUC,FullRocBlindAUC);
					if (rocadded == 0)
					{
						split.blindSen <- splitRoc;
					}
					else
					{
						split.blindSen <- rbind(split.blindSen,splitRoc);
					}
					rocadded = rocadded + 1;
				}
				
				totsize <- totsize + framesize;
				scase <- sum(BlindSet[,Outcome] == 1);
				scotr <- sum(BlindSet[,Outcome] == 0);
				sizecases <- sizecases + scase;
				sizecontrol <- sizecontrol + scotr;
				psen <- sum( 1*((BlindSet[,Outcome] > 0)*( p >= 0.0 )) , na.rm = TRUE)
				pspe <- sum( 1*((BlindSet[,Outcome] == 0)*( p < 0.0 )) , na.rm = TRUE)
				acc <- acc + psen + pspe;
				sen <- sen + psen;
				spe <- spe + pspe;
				psen <- sum( 1*((BlindSet[,Outcome] > 0)*( Full.p >= 0.0 )) , na.rm = TRUE)
				pspe <- sum( 1*((BlindSet[,Outcome] == 0)*( Full.p < 0.0 )) , na.rm = TRUE)
				Full.acc <- Full.acc + psen + pspe;
				Full.sen <- Full.sen + psen;
				Full.spe <- Full.spe + pspe;
				paracc = acc/totsize;
				psen = 0;
				pspe = 0;
				if (sizecases>0) 
				{
					psen = sen/sizecases;
				}
				if (sizecontrol>0) 
				{
					pspe = spe/sizecontrol;
				}

				Full.paracc = Full.acc/totsize;
				Full.psen = 0;
				Full.pspe = 0;
				if (sizecases>0) 
				{
					Full.psen = Full.sen/sizecases;
				}
				if (sizecontrol>0) 
				{
					Full.pspe = Full.spe/sizecontrol;
				}

				WholeFoldBlindAccuracy <- append(WholeFoldBlindAccuracy,redBootCross$blind.accuracy);
				WholeFoldBlindSpecificity <- append(WholeFoldBlindSpecificity,redBootCross$blind.specificity);
				WholeFoldBlindSensitivity <- append(WholeFoldBlindSensitivity,redBootCross$blind.sensitivity);

				FoldBlindAccuracy <- append(FoldBlindAccuracy,redBootCross_S$blind.accuracy);
				FoldBlindSpecificity <- append(FoldBlindSpecificity,redBootCross_S$blind.specificty);
				FoldBlindSensitivity <- append(FoldBlindSensitivity,redBootCross_S$blind.sensitivity);

				Full.ptrain <- predictForFresa(redBootCross$boot.model,TrainSet, 'linear');
				ptrain <- predictForFresa(redfoldmodel,TrainSet, 'linear');

				if ( cycleinsert == 1)
				{
					cvcycle.predictions <- cbind(BlindSet[,Outcome],p.BBH,i);
				}
				if (inserted == 1)
				{
					totSamples <- cbind(BlindSet[,Outcome],p,i,medianPred,baggedForwardPredict,p.forward,p.BBH);
					rownames(totSamples) <- rownames(BlindSet);
					Full.totSamples <- cbind(BlindSet[,Outcome],Full.p,i);
					rownames(Full.totSamples) <- rownames(BlindSet);
					KNNSamples <- cbind(BlindSet[,Outcome],abs(knnclass$prob$prob-1*(knnclass$prediction=="0")),i);
					rownames(KNNSamples) <- rownames(BlindSet);
					Full.KNNSamples <- cbind(BlindSet[,Outcome],abs(Fullknnclass$prob$prob-1*(Fullknnclass$prediction=="0")),i);
					rownames(Full.KNNSamples) <- rownames(BlindSet);
					totTrainSamples <- cbind(TrainSet[,Outcome],ptrain,i);
					rownames(totTrainSamples) <- rownames(TrainSet);
					Full.totTrainSamples <- cbind(TrainSet[,Outcome],Full.ptrain,i);
					rownames(Full.totTrainSamples) <- rownames(TrainSet);
					
				}
				else
				{
					px <- cbind(BlindSet[,Outcome],p,i,medianPred,baggedForwardPredict,p.forward,p.BBH);
					rownames(px) <- rownames(BlindSet);
					totSamples <- rbind(totSamples,px);
					px <- cbind(BlindSet[,Outcome],p.BBH,i);
					rownames(px) <- rownames(BlindSet);
					cvcycle.predictions <- rbind(cvcycle.predictions,px);
					px <- cbind(BlindSet[,Outcome],Full.p,i);
					rownames(px) <- rownames(BlindSet);
					Full.totSamples <- rbind(Full.totSamples,px);
					px <- cbind(BlindSet[,Outcome],abs(knnclass$prob$prob-1*(knnclass$prediction=="0")),i);
					rownames(px) <- rownames(BlindSet);
					KNNSamples <- rbind(KNNSamples,px);
					px <- cbind(BlindSet[,Outcome],abs(Fullknnclass$prob$prob-1*(Fullknnclass$prediction=="0")),i);
					rownames(px) <- rownames(BlindSet);
					Full.KNNSamples <- rbind(Full.KNNSamples,px);
					px <- cbind(TrainSet[,Outcome],ptrain,i);
					rownames(px) <- rownames(TrainSet);
					totTrainSamples <- rbind(totTrainSamples,px);
					px <- cbind(TrainSet[,Outcome],Full.ptrain,i);
					rownames(px) <- rownames(TrainSet);
					Full.totTrainSamples <- rbind(Full.totTrainSamples,px);
				}
				
				
				formulas <- append(formulas,as.character(redCurmodel_S$back.formula)[3]);
				BeforeBHFormulas <- append(BeforeBHFormulas,as.character(redCurmodel_S$beforeFSC.formula)[3]);
				ForwardFormulas <- append(ForwardFormulas,as.character(UCurModel_S$formula)[3]);
				baggFormulas <- append(baggFormulas,bagg$formula);

				knnACC <- sum(KNNSamples[,1] == (KNNSamples[,2]>0.5))/totsize;
				knnSEN <- sum((KNNSamples[,1]>0.5) & (KNNSamples[,2]>0.5))/sizecases;
				knnSPE <- sum((KNNSamples[,1]<0.5) & (KNNSamples[,2]<0.5))/sizecontrol;

				Full.knnACC <- sum(Full.KNNSamples[,1] == (Full.KNNSamples[,2]>0.5))/totsize;
				Full.knnSEN <- sum((Full.KNNSamples[,1]>0.5) & (Full.KNNSamples[,2]>0.5))/sizecases;
				Full.knnSPE <- sum((Full.KNNSamples[,1]<0.5) & (Full.KNNSamples[,2]<0.5))/sizecontrol;


				cat ("Loop :",i,"Blind Cases =",scase,"Blind Control =",scotr,"Total =",totsize, "Size Cases =",sizecases,"Size Control =",sizecontrol,"\n")
				cat ("Accumulated Models CV Accuracy        =",paracc,"Sensitivity =",psen,"Specificity =",pspe,"\n")
				cat ("Initial Model Accumulated CV Accuracy =",Full.paracc,"Sensitivity =",Full.psen,"Specificity =",Full.pspe,"\n");
				cat ("Initial Model Bootstrapped Accuracy   =",redBootCross$blind.accuracy,"Sensitivity =",redBootCross$blind.sensitivity,"Specificity =",redBootCross$blind.specificity,"\n")
				cat ("Current Model Bootstrapped Accuracy   =",redBootCross_S$blind.accuracy,"Sensitivity =",redBootCross_S$blind.sensitivity,"Specificity =",redBootCross_S$blind.specificity,"\n")
				cat ("Current KNN Accuracy   =",knnACC,"Sensitivity =",knnSEN,"Specificity =",knnSPE,"\n")
				cat ("Initial KNN Accuracy   =",Full.knnACC,"Sensitivity =",Full.knnSEN,"Specificity =",Full.knnSPE,"\n")
				cat ("Train Correlation: ",tcor," Blind Correlation :",bcor,"\n KNN to Model Confusion Matrix: \n")
				print(table(KNNSamples[,2]>0.5,totSamples[,2]>0.0))
			}
			# else
			# {
				# cat ("Loop :",i,"No Model.\n")
			# }
			
		}
		else
		{
			cat ("Loop :",i,"No Model.\n")
		}

		uniEval <- getVar.Bin(Full_formula,TrainSet,Outcome,type = type,testData=BlindSet);
		if (i==1)
		{
			uniTrainAccuracy <- rbind(uniEval$uniTrainAccuracy);
			TopUniTrainCor <- vector();
		}
		else
		{
			uniTrainAccuracy <- rbind(uniTrainAccuracy,uniEval$uniTrainAccuracy);
		}
		if ( j == 1)
		{
			cvcycle.uniAccuracies <- uniEval$uniTestAccuracy * framesize;
			totblindadded = framesize;
			topUniTestCor <- vector();
			totalUniCor = 0; 
		}
		else
		{
			cvcycle.uniAccuracies <- rbind(cvcycle.uniAccuracies,uniEval$uniTestAccuracy * framesize);
			totblindadded = totblindadded + framesize;
		}

		if (lastinserted<inserted)
		{		
			uniEvalCor <- getVar.Bin(redCurmodel_S$back.formula,TrainSet,Outcome,type = type,testData=BlindSet);
			TopUniTrainCor <- append(TopUniTrainCor,uniEvalCor$uniTrainAccuracy[1]);
			topUniTestCor <- append(topUniTestCor,uniEvalCor$uniTestAccuracy[1] * framesize);
			totalUniCor <- totalUniCor + framesize
		}
		
				
		if ( j == K)
		{
			if (totalUniCor>0) TopUniCoherenceTest <- append(TopUniCoherenceTest,sum(topUniTestCor)/totalUniCor)
			if (i == K)
			{
				uniTestAccuracy <- rbind(colSums(cvcycle.uniAccuracies)/totblindadded);
			}
			else
			{
				uniTestAccuracy <- rbind(uniTestAccuracy,colSums(cvcycle.uniAccuracies)/totblindadded);
			}
		}



		if ( j == K)
		{
			nsamp <- nrow(cvcycle.predictions)
			if (nsamp>0)
			{
				atRocAUC <- pROC::roc(cvcycle.predictions[,1], cvcycle.predictions[,2],plot=FALSE,auc=TRUE,smooth=FALSE)$auc;
				testAccuracy <- append(testAccuracy,sum(cvcycle.predictions[,1] == 1.0*(cvcycle.predictions[,2]>=0.0))/nsamp);
				testSensitivity <- append(testSensitivity,sum((cvcycle.predictions[,1] == 1) & (cvcycle.predictions[,2]>=0.0))/sum(cvcycle.predictions[,1] == 1));
				testSpecificity <- append(testSpecificity,sum((cvcycle.predictions[,1] == 0) & (cvcycle.predictions[,2] <0.0))/sum(cvcycle.predictions[,1] == 0));
				testAUC <- append(testAUC,atRocAUC);
			}
#			print(testAccuracy)
#			print(testAUC)
		}
		
	}
	if (length(formulas)==0)
	{
		stop("No Significant Models Found\n");
	}
	colnames(totSamples) <- c("Outcome","Prediction","Model","Median","Bagged","Forward","Backwards");
	totSamples <- as.data.frame(totSamples);

	colnames(Full.totSamples) <- c("Outcome","Prediction","Model");
	Full.totSamples <- as.data.frame(Full.totSamples);

	colnames(totTrainSamples) <- c("Outcome","Prediction","Model");
	totTrainSamples <- as.data.frame(totTrainSamples);
	colnames(Full.totTrainSamples) <- c("Outcome","Prediction","Model");
	Full.totTrainSamples <- as.data.frame(Full.totTrainSamples);

	colnames(KNNSamples) <- c("Outcome","Prediction","Model");
	KNNSamples <- as.data.frame(KNNSamples);
	
	colnames(Full.KNNSamples) <- c("Outcome","Prediction","Model");
	Full.KNNSamples <- as.data.frame(Full.KNNSamples);

	
	
	BSWiMS.ensemble.prediction <- NULL

	bsta <- boxplot(totSamples$Prediction~rownames(totSamples),plot=FALSE)
	sta <- cbind(bsta$stats[3,])
	rownames(sta) <- bsta$names
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

	
	
	plotModels.ROC(totSamples);
	par(mfrow=c(1,1))
	incBsen=0
	aucBlindTest <- pROC::roc(totSamples[,1],totSamples[,2],col="red",auc=TRUE,plot=TRUE,smooth=FALSE,lty=3)$auc
	par(new=TRUE)
	aucCVBlind <- pROC::roc(Full.totSamples[,1],Full.totSamples[,2],col="blue",auc=TRUE,plot=TRUE,ci=TRUE,smooth=FALSE)$auc
	par(new=TRUE)
	aucTrain <- pROC::roc( FullBootCross$outcome, FullBootCross$boot.model$linear.predictors,col="green",plot=TRUE,auc=TRUE,smooth=FALSE)$auc;        
	par(new=TRUE)
	aucBoot <- pROC::roc( FullBootCross$testOutcome, FullBootCross$testPrediction,col="black",auc=TRUE,plot=TRUE,smooth=FALSE)$auc;
	ley.names <- c(paste("Bootstrapped: Train Model ROC (",sprintf("%.3f",aucTrain),")"),paste("Bootstrapped: Blind ROC (",sprintf("%.3f",aucBoot),")"),
	paste("CV: Blind ROC (",sprintf("%.3f",aucCVBlind),")"),paste("CV: Blind Fold Models Coherence (",sprintf("%.3f",aucBlindTest),")"))
	ley.colors <- c("green","black","blue","red")
	ley.lty <- c(1,1,1,3)
	if (rocadded>0)
	{
		boxplot(split.blindSen,add=TRUE, axes = FALSE,boxwex=0.04,at=specificities);
		sumSen <- colMeans(split.blindSen,na.rm = TRUE);
		sennames <- names(sumSen);
		sumSen <- append(0,sumSen);
		sumSen <- append(sumSen,1);
		sennames <- append("1",sennames);
		sennames <- append(sennames,"0");
		names(sumSen) <- sennames;
		spevalues <- as.numeric(names(sumSen));
		lines(spevalues,sumSen,col="red",lwd=2.0);
		auc = 0;
		for (i in 2:length(spevalues))
		{
			auc = auc + (spevalues[i-1]-spevalues[i])*(sumSen[i-1]+(sumSen[i]-sumSen[i-1])/2)
		}
		ley.names <- append(ley.names,paste("CV Blind: Mean ROC of Models (",sprintf("%.3f",auc),")"));
		ley.colors <- append(ley.colors,"red");
		ley.lty  <- append(ley.lty,1);
	}
	else
	{
		sumSen = NA;
	}
	
	legend(0.6,0.30, legend=ley.names,col = ley.colors, lty = ley.lty,bty="n")


	if (!is.null(uniTrainAccuracy))
	{
		uniTrainAccuracy <- as.data.frame(uniTrainAccuracy);
		uniTestAccuracy <- as.data.frame(uniTestAccuracy);
		colnames(uniTrainAccuracy) <-  attr(terms(redCurmodel_Full$back.formula),'term.labels');
		colnames(uniTestAccuracy) <-  attr(terms(redCurmodel_Full$back.formula),'term.labels');
	}
	
	result <- list(formula.list=formulas,
	Models.testPrediction=totSamples,
	FullBWiMS.testPrediction=Full.totSamples,
	TestRetrained.blindPredictions=blindreboot,
	LastTrainBSWiMS.bootstrapped=redCurmodel_S$bootCV,
	Test.accuracy=paracc,
	Test.sensitivity=psen,
	Test.specificity=pspe,
	Train.correlationsToFull=trainCorrelations,
	Blind.correlationsToFull=blindCorrelations,
	FullModelAtFoldAccuracies=WholeFoldBlindAccuracy,
	FullModelAtFoldSpecificties=WholeFoldBlindSpecificity,
	FullModelAtFoldSensitivities=WholeFoldBlindSensitivity,
	FullModelAtFoldAUC=WholeFoldBlindAUC,
	CVTrain.Accuracies=trainAccuracy,
	CVTrain.Sensitivity=trainSensitivity,
	CVTrain.Specificity=trainSpecificity,
	CVTrain.AUCs=trainAUC,
	CVTest.Accuracies=testAccuracy,
	CVTest.Sensitivity=testSensitivity,
	CVTest.Specificity=testSpecificity,
	CVTest.AUCs=testAUC,
	AtCVFoldModelBlindAccuracies=FoldBlindAccuracy,
	AtCVFoldModelBlindSpecificities=FoldBlindSpecificity,
	AtCVFoldModelBlindSensitivities=FoldBlindSensitivity,
	Models.CVblindMeanSensitivites=sumSen,
	forwardSelection = CurModel_Full,
	updateforwardSelection = UCurModel_Full,
	BiSWiMS = redCurmodel_Full,
	FullBWiMS.bootstrapped=FullBootCross,
	Models.testSensitivities = split.blindSen,
	FullKNN.testPrediction=Full.KNNSamples,
	KNN.testPrediction=KNNSamples,
	Fullenet=Fullenet,
	LASSO.testPredictions=enetSamples,
	LASSOVariables=LASSOVariables,
	uniTrain.Accuracies=uniTrainAccuracy,
	uniTest.Accuracies=uniTestAccuracy,
	uniTest.TopCoherence=TopUniCoherenceTest,
	uniTrain.TopCoherence=TopUniTrainCor,
	Models.trainPrediction=totTrainSamples,
	FullBWiMS.trainPrediction=Full.totTrainSamples,
	LASSO.trainPredictions=enetTrainSamples,
	BSWiMS.ensemble.prediction = BSWiMS.ensemble.prediction,
	ForwardFormulas.list = ForwardFormulas,
	BeforeBHFormulas.list = BeforeBHFormulas,
	baggFormulas.list = baggFormulas
	);
	return (result)
}
