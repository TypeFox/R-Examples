bootstrapValidation_Bin <-
function(fraction=1.00,loops=200,model.formula,Outcome,data,type=c("LM","LOGIT","COX"),plots=TRUE)
{

	if (class(model.formula) != "formula")
	{
		model.formula = formula(model.formula);
	}

	varsList <- as.list(attr(terms(model.formula),"variables"))
	if (length(varsList)>2)
	{

		 type <- match.arg(type);
		 casesample = subset(data,get(Outcome)== 1);
		 controlsample = subset(data,get(Outcome) == 0);
		 sizecases = nrow(casesample);
		 sizecontrol = nrow(controlsample);

		modelMat <- model.matrix(model.formula,data);
		if (type=="COX")
		{
			response <- data[,all.vars(model.formula)[1:2]];
		}
		else
		{
			response <- data[,Outcome];
		}

		output<-.Call("bootstrapValidationBinCpp", fraction, loops, modelMat,type, data.matrix(response));

		basemodel <- modelFitting(model.formula,data,type)
		colnames(output$bcoef) <- names(basemodel$coefficients);
		basepredict <- predictForFresa(basemodel,testData=data,predictType = 'linear');
		basemodel$linear.predictors <- basepredict;
		framesize = nrow(data);
		acc = 0.0;
		sen = 0.0;
		spe = 0.0;
		for (i in 1:framesize)
		{
			if ((data [i,Outcome] > 0) && (basepredict[i] > 0) ) 
			{
				acc = acc + 1.0; 
				sen = sen + 1.0;
			}
			if ((data [i,Outcome] == 0) && (basepredict[i] < 0) ) 
			{
				acc = acc + 1.0; 
				spe = spe + 1.0;
			}
		}
		baseAcc = acc/framesize;
		baseSen = sen/sizecases;
		baseSpe = spe/sizecontrol;
		bootmodel <- basemodel;
		wt <- 2.0*(output$accuracy-0.5); 			# the linear weights
		wt <- 0.1+wt*(wt > 0);

		bootmodel$coefficients <-  apply(output$bcoef,2,weighted.mean,wt = wt, na.rm = TRUE)


		if (length(output$sumwts)>0)
		{
			for (n in 1:length(output$sumwts))
			{
				if (!is.na(output$sumwts[n]))
				{
					if (output$sumwts[n] < 0)  
					{
						output$sumwts[n] = 1;
						output$sumwtdcf[n] = 0;  
					}
				}
				else
				{
					output$sumwts[n] = 1;
					output$sumwtdcf[n] = 0;  
				}
			}
		}

		pr <- predictForFresa(bootmodel,testData=data,predictType = 'linear');
		bootmodel$linear.predictors <- pr;
		p <- predictForFresa(bootmodel,testData=data,predictType = 'prob');
		bootmodel$fitted.values <- p;
		sen <- sum( 1*((data[,Outcome] > 0)*( bootmodel$linear.predictors >= 0 )) , na.rm = TRUE)
		spe <- sum( 1*((data[,Outcome] == 0)*( bootmodel$linear.predictors < 0 )) , na.rm = TRUE)
		acc <- sen+spe;
		acc = acc/framesize;
		sen = sen/sizecases;
		spe = spe/sizecontrol;
		blidRoc = NULL;
		bootRoc = NULL;
		if (plots)
		{
			par(mfrow=c(2,2))
			pROC::roc( data[,Outcome], basemodel$linear.predictors, col="blue",plot=TRUE,smooth=FALSE,progress= 'none');
			par(new=TRUE)
			blidRoc <- pROC::roc(output$testoutcome,output$testprediction,col="red",auc=TRUE,print.auc=TRUE,plot=TRUE,smooth=FALSE,progress= 'none')
			par(new=TRUE)
			bootRoc <- pROC::roc( data[,Outcome], bootmodel$linear.predictors,plot=TRUE,ci=plots,auc=TRUE,of='se',specificities=c(0.95,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.05),boot.n=200,smooth=FALSE,progress= 'none');
				plot(ecdf(output$taccuracy),main="Accuracy",xlim=c(0.5, 1.0),col="black",lwd = 2,verticals = TRUE, do.points = FALSE);
				abline(v=output$BlindAccuracy,col = "red");
				abline(v=acc,col = "blue");
				plot(ecdf(output$tsensitivity),main="Sensitivity",xlim=c(0.5, 1.0),col="black",lwd = 2,verticals = TRUE, do.points = FALSE);
				abline(v=output$BlindSensitivity,col = "red");
				abline(v=sen,col = "blue");
				plot(ecdf(output$tspecificity),main="Specificity",xlim=c(0.5, 1.0),col="black",lwd = 2,verticals = TRUE, do.points = FALSE);
				abline(v=output$BlindSpecificity,col = "red");
				abline(v=spe,col = "blue");
			par(mfrow=c(1,1))
#			cat("Mean Wt:",mean(wt)," Max Wt:",max(wt)," Min wt:",min(wt),"\n")
		}
		else
		{
			blidRoc <- llist(auc=output$resul$blindAUC);
			bootRoc <- pROC::roc( data[,Outcome], bootmodel$linear.predictors,plot=FALSE,auc=TRUE,smooth=FALSE,progress= 'none');
		}
		colnames(output$IDI) <- attr(terms(model.formula),'term.labels');
		colnames(output$zIDI) <- attr(terms(model.formula),'term.labels');
		colnames(output$NRI) <- attr(terms(model.formula),'term.labels');
		colnames(output$zNRI) <- attr(terms(model.formula),'term.labels');
		colnames(output$bcoef) <- names(basemodel$coefficients);
		tem=output$resul;

		result <- structure(llist(
		data=data,
		outcome=data[,Outcome],
		blind.accuracy=output$BlindAccuracy,
		blind.sensitivity=output$BlindSensitivity,
		blind.specificity=output$BlindSpecificity,
		train.ROCAUC=output$trainRoc,
		blind.ROCAUC= blidRoc,
		boot.ROCAUC=bootRoc, 
		fraction=fraction,
		loops=loops,
		base.Accuracy=baseAcc,
		base.sensitivity=baseSen,
		base.specificity=baseSpe,
		accuracy=output$accuracy,
		sensitivity=output$sensitivity,
		specificity=output$specificity,
		train.accuracy=output$taccuracy,
		train.sensitivity=output$tsensitivity,
		train.specificity=output$tspecificity,
		s.coef=output$bcoef,
		boot.model=bootmodel,
		boot.accuracy=acc,
		boot.sensitivity=sen,
		boot.specificity=spe,
		z.NRIs=output$zNRI,
		z.IDIs=output$zIDI,
		test.z.NRIs=tem$test_zNRI,
		test.z.IDIs=tem$test_zIDI,
		startBetas=tem$rawbetas,
		NRIs=output$NRI,
		IDIs=output$IDI,
		testOutcome=output$testoutcome,
		testPrediction=output$testprediction,
		labels = FALSE), class = "bootstrapValidation_Bin");
	}
	else
	{
		result <- NULL;
	}
	return (result);
}

