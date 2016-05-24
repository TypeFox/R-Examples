getVar.Bin <-
function (object,data,Outcome="Class", type = c("LOGIT", "LM","COX"),testData=NULL,callCpp=TRUE) 
{

	if (is.null(testData))
	{
		testData <- data;
	}

	model_zidi <- vector();
	model_znri <- vector();
	model_idi <- vector();
	model_nri <- vector();

	t.model_zidi <- vector();
	t.model_znri <- vector();
	t.model_idi <- vector();
	t.model_nri <- vector();
	test.accuracy <- vector();
	train.accuracy <- vector();
  
	varsList <- as.list(attr(terms(object),"variables"))

		
	outCome = paste(varsList[2]," ~ 1");
	startlist = 3;
	frm1 = outCome;
	for ( i in startlist:length(varsList))
	{
		frm1 <- paste(frm1,paste(" + ",varsList[i]));
	}
	ftmp <- formula(frm1);
	sizetest <- nrow(testData);
	sizetrain <- nrow(data);
	
	FullModel <- modelFitting(ftmp,data,type,fast=callCpp)
	if ( !inherits(FullModel, "try-error"))
	{
	
		FullPredict_train <- predictForFresa(FullModel,data,'prob');
		FullPredict <- predictForFresa(FullModel,testData, 'prob');

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
			if (inherits(redModel, "try-error"))
			{
				redModel <- FullModel;
			}

			redPredict <- predictForFresa(redModel,testData,'prob');
			redPredict_train <- predictForFresa(redModel,data,'prob');
			
			if (!callCpp)
			{
				iprob <- improveProb(redPredict,FullPredict,testData[,Outcome]);
				iprob_t <- improveProb(redPredict_train,FullPredict_train,data[,Outcome]);
			}
			else
			{
				iprob <-.Call("improveProbCpp",redPredict,FullPredict,testData[,Outcome],sizetrain);
				iprob_t <-.Call("improveProbCpp",redPredict_train,FullPredict_train,data[,Outcome],0);
			}

#univariate analysis
			uniModel <- modelFitting(formula(paste(outCome," + ",varsList[i])),data,type,fast=callCpp);
			uniPredict <- predictForFresa(uniModel,testData,'prob');
			uniPredict_train <- predictForFresa(uniModel,data,'prob');
			test.accuracy = append(test.accuracy,sum(testData[,Outcome] == 1.0*(uniPredict>=0.5))/sizetest);
			train.accuracy = append(train.accuracy,sum(data[,Outcome] == 1.0*(uniPredict_train>=0.5))/sizetrain);
			
#end univariate analysis
			
			model_zidi <- append(model_zidi,iprob$z.idi);
			model_idi <- append(model_idi,iprob$idi);
			model_nri <- append(model_nri,iprob$nri);
			model_znri <- append(model_znri,iprob$z.nri);

			t.model_zidi <- append(t.model_zidi,iprob_t$z.idi);
			t.model_idi <- append(t.model_idi,iprob_t$idi);
			t.model_nri <- append(t.model_nri,iprob_t$nri);
			t.model_znri <- append(t.model_znri,iprob_t$z.nri);
		}

	}
	
	 result <- list(
		 z.IDIs=t.model_zidi,
		 z.NRIs=t.model_znri,
		 IDIs=t.model_idi,
		 NRIs=t.model_nri,
		 testData.z.IDIs=model_zidi,
		 testData.z.NRIs=model_znri,
		 testData.IDIs=model_idi,
		 testData.NRIs=model_nri,
		 uniTrainAccuracy = train.accuracy,
		 uniTestAccuracy = test.accuracy
	 );

    return (result)
}
