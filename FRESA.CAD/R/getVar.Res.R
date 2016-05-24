getVar.Res <-
function (object,data,Outcome="Class", type=c("LM","LOGIT","COX"),testData=NULL,callCpp=TRUE) 
{
	if (is.null(testData))
	{
		testData <- data;
	}
  
	samples <- nrow(data);
	varsList <- as.list(attr(terms(object),"variables"))
	termList <- as.list(attr(terms(object),"term.labels"))
	
	outCome = paste(varsList[2]," ~ 1");
	startlist = 3;
	frm1 = outCome;
	for ( i in startlist:length(varsList))
	{
		frm1 <- paste(frm1,paste(" + ",varsList[i]));
	}
	ftmp <- formula(frm1);
	FullModel <- modelFitting(ftmp,data,type,fast=callCpp);

	FullResiduals <- residualForFRESA(FullModel,data,Outcome);
	testResiduals <- residualForFRESA(FullModel,testData,Outcome);

	model_tpvalue <- vector();
	model_bpvalue <- vector();
	model_neri <- vector();
	model_wpvalue <- vector();
	model_fpvalue <- vector();
	testmodel_tpvalue <- vector();
	testmodel_bpvalue <- vector();
	testmodel_neri <- vector();
	testmodel_wpvalue <- vector();
	testmodel_fpvalue <- vector();
	unitest.MSS <- vector();
	unitrain.MSS <- vector();
	if (length(termList)>0)
	{
		for ( i in startlist:length(varsList))
		{
		
			frm1 = outCome;
			for ( j in startlist:length(varsList))
			{
				if (i!=j)
				{
					frm1 <- paste(frm1,paste(" + ",varsList[j]));
				}
			}
			ftmp <- formula(frm1);
			redModel <- modelFitting(ftmp,data,type,fast=callCpp)

			if ( inherits(redModel, "try-error"))
			{
				redModel <- FullModel;
			}
			
			redResiduals <- residualForFRESA(redModel,data,Outcome);
			redTestResiduals <- residualForFRESA(redModel,testData,Outcome);

			if (nrow(testData)>2) # at least 3 samples for improvement analysis
			{
				if (!callCpp)
				{
					iprob <- improvedResiduals(redResiduals,FullResiduals);
					testiprob <- improvedResiduals(redTestResiduals,testResiduals);
				}
				else
				{
					iprob <- .Call("improvedResidualsCpp",redResiduals,FullResiduals," ",0);
					testiprob <- .Call("improvedResidualsCpp",redTestResiduals,testResiduals," ",samples);
				}
				model_tpvalue <- append(model_tpvalue,iprob$tP.value);
				model_bpvalue <- append(model_bpvalue,iprob$BinP.value);
				model_wpvalue <- append(model_wpvalue,iprob$WilcoxP.value);
				model_fpvalue <- append(model_fpvalue,iprob$FP.value);
				model_neri <- append(model_neri,iprob$NeRI);
				testmodel_tpvalue <- append(testmodel_tpvalue,testiprob$tP.value);
				testmodel_bpvalue <- append(testmodel_bpvalue,testiprob$BinP.value);
				testmodel_wpvalue <- append(testmodel_wpvalue,testiprob$WilcoxP.value);
				testmodel_fpvalue <- append(testmodel_fpvalue,testiprob$FP.value);
				testmodel_neri <- append(testmodel_neri,testiprob$NeRI);
			}

#univariate analysis
			uniModel <- modelFitting(formula(paste(outCome," + ",varsList[i])),data,type,fast=callCpp);
			uniPredict <- predictForFresa(uniModel,testData,'linear');
			uniPredict_train <- predictForFresa(uniModel,data,'linear');
			unitest.MSS = append(unitest.MSS,mean((uniPredict-testData[,Outcome])^2));
			unitrain.MSS = append(unitrain.MSS,mean((uniPredict_train-data[,Outcome])^2));
#end univariate analysis
			
		}
	}
	
	 result <- list(
	 tP.value=model_tpvalue,
	 BinP.value=model_bpvalue,
	 WilcoxP.value=model_wpvalue,
	 FP.value=model_fpvalue,
	 NeRIs=model_neri,
	 testData.tP.value=testmodel_tpvalue,
	 testData.BinP.value=testmodel_bpvalue,
	 testData.WilcoxP.value=testmodel_wpvalue,
	 testData.FP.value=testmodel_fpvalue,
	 testData.NeRIs=testmodel_neri,
	 unitestMSS = unitest.MSS,
	 unitrainMSS = unitrain.MSS
	 );

    return (result)
}
