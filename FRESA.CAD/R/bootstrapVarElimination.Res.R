bootstrapVarElimination_Res <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),loops=250,fraction=1.00,setIntersect=1,print=TRUE,plots=TRUE,adjsize=1) 
{
  	testType <- match.arg(testType)

#	boot.var.NeRISelection <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),loops,fraction=1.0,setIntersect=1,NeRICV=NULL) 
	boot.var.NeRISelection <- function (objectt,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),loops,fraction=1.0,setIntersect=1) 
	{
		testType <- match.arg(testType)
		type <- match.arg(type);
		FullModel <- objectt;
		varsList <- as.list(attr(terms(objectt),"variables"))
		
		climpvalue = 0.4999; # 
		removeID = 0;

		outCome = paste(varsList[2]," ~ ",setIntersect);
		startlist = 3 ;
		frm1 = outCome;
		for ( i in startlist:length(varsList))
		{
			frm1 <- paste(frm1,paste(" + ",varsList[i]));
		}
		ftmp <- formula(frm1);
		backfrm <- frm1;
#		cat("Start  Formula :",frm1,"\n")
		NeRICVp <- bootstrapValidation_Res(fraction,loops,ftmp,Outcome,data,type,plots=plots)
#		cat("Bootstrapped  Formula :",frm1,"\n")
		startSearch = startlist + startOffset;
		maxPvalue=pvalue;
		if (length(varsList) >= startSearch)
		{		
			idlist=startOffset+1;
			frm1 = outCome;
			if ((startSearch-1) >= startlist)
			{
				for ( i in startlist:(startSearch-1))
				{
					frm1 <- paste(frm1,paste(" + ",varsList[i]));
				}
			}
#			cat("Elimination Base Formula :",frm1,"\n")
			idx = 2;
			who = -1;
			maxPvalue = pvalue;
			for ( i in startSearch:length(varsList))
			{
				if (any(is.na(NeRICVp$bin.pvlaues[,idlist])))
				{
					who = i;
				}
				else
				{ # reduce probability by two to test for equivalence of reduced model to the Full model
					switch(testType, 
						tStudent = 
						{ 
							ci <- as.vector(quantile(NeRICVp$tStudent.pvalues[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
							ci2 <- as.vector(quantile(NeRICVp$test.tStudent.pvalues[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
						},
						Wilcox = 
						{ 
							ci <- as.vector(quantile(NeRICVp$wilcox.pvalues[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
							ci2 <- as.vector(quantile(NeRICVp$test.wilcox.pvalues[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
						},
						Binomial =
						{ 
							ci <- as.vector(quantile(NeRICVp$bin.pvlaues[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
							ci2 <- as.vector(quantile(NeRICVp$test.bin.pvlaues[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
						},
						Ftest =
						{ 
							ci <- as.vector(quantile(NeRICVp$F.pvlaues[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
							ci2 <- as.vector(quantile(NeRICVp$test.F.pvlaues[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
						},
					)
					cmax = max(ci[2],ci2[2]);
#					cmax = ci2[2];
					if (cmax >= maxPvalue)
					{
						maxPvalue = cmax;
						who = i;
					}
					if (ci[3] > maxPvalue)
					{
						who = i;
					}
				}
				idlist=idlist+1;
			}
			for ( i in startSearch:length(varsList))
			{
				if (who != i)
				{
					frm1 <- paste(frm1,paste(" + ",varsList[i]));
				}
				else
				{
					removeID=idlist;
#					cat ("Removed ",paste(" -> ",varsList[i]),"\n");
				}
			}
#			cat ("Rm: Formula: ",frm1,"\n")
			if ((length(varsList) == startSearch) && (who == startSearch)) 
			{
				removeID = -1;
			}
			ftmp <- formula(frm1);
			FullModel <- modelFitting(ftmp,data,type)
			backfrm <- frm1
			if (inherits(FullModel, "try-error"))
			{
#				cat("Error: Reduced Formula: ",frm1,"\n");
				who= length(varsList)
				frm1 = outCome;
				if ((startSearch-1) >= startlist)
				{
					for ( i in startlist:(startSearch-1))
					{
						frm1 <- paste(frm1,paste(" + ",varsList[i]));
					}
				}
				for ( i in startSearch:length(varsList))
				{
					if (who != i)
					{
						frm1 <- paste(frm1,paste(" + ",varsList[i]));
					}
				}
				ftmp <- formula(frm1);
				backfrm <- frm1
				FullModel <- modelFitting(ftmp,data,type);
			}

		}
		NeRICVp <- bootstrapValidation_Res(fraction,loops,ftmp,Outcome,data,type,plots=plots)
		testRMSE <- NeRICVp$testRMSE;
		if (length(testRMSE)==0) testRMSE=1.0e10;
		if (is.na(testRMSE)) testRMSE=1.0e10;
		if (is.null(testRMSE)) testRMSE=1.0e10;
		result <- list(Removed=removeID,backfrm=backfrm,testRMSE=testRMSE,max.pvalue=maxPvalue);
		
		return (result)
	}

	model <- NULL;
	NeRICV <- NULL;
	modelReclas <- NULL;
	minbootRMSE <- 1.0;
	stopFratio = 1.025; # stop if there is a small increase in RMS test error 
	changes=1;
	if (adjsize>1)
	{
		bkobjt <- bootstrapVarElimination_Res(object,pvalue,Outcome,data,startOffset,type,testType,loops,fraction,setIntersect,print,plots,adjsize=1); #just remove the features that produce similar models
#		object <- bkobjt$bootCV$boot.model; 
		object <- modelFitting(bkobjt$back.formula,data,type);
		adjsize = floor(adjsize);	
#		cat("Adjust size:",adjsize,"\n");
		stopFratio = 1.10; # max increase in RMS test error	
		if (is.null(bkobjt$bootCV)) changes=0;
	}
	modelReclas <-getVar.Res(object,data=data,Outcome=Outcome,type);
	NeRICV <- bootstrapValidation_Res(fraction,loops,formula(object),Outcome,data,type,plots=plots)

	minbootRMSE <- NeRICV$testRMSE;		
	startRMSE <- minbootRMSE;

	varsList <- as.list(attr(terms(object),"variables"))
	outCome = paste(varsList[2]," ~ ",setIntersect);
	startlist = 3 ;
	frm1 = outCome;
	if (startlist <= length(varsList))
	{
		for ( i in startlist:length(varsList))
		{
			frm1 <- paste(frm1,paste(" + ",varsList[i]));
		}
	}
	model.formula <- frm1;
	
	
	loopsAux=0;
    model = object;
	beforeFSCmodel <- object;
	beforeFSC.model.formula <- frm1;
	wts <- rep(1,length(beforeFSCmodel$coefficients));
	names(wts) <- names(beforeFSCmodel$coefficients);
	changes2 <- 0
	changed <- 0;
	bk <- NULL;
	while ((changes>0) && (loopsAux<100)) 
	{
		p.elimin <- pvalue;
		if (adjsize>1)
		{
			modsize <- length(as.list(attr(terms(model),"term.labels")));
			if (modsize<1) modsize=1;
			qvalue <- 4*pvalue;
			if (qvalue < 0.1) qvalue=0.1; # lets keep the minimum q value to 0.1
			p.elimin <- min(pvalue,modsize*qvalue/adjsize) # # BH alpha the elimination p-value
		}

		bk = boot.var.NeRISelection(objectt=model,pvalue=p.elimin,Outcome=Outcome,data=data,
		startOffset=startOffset,type=type,testType=testType,loops=loops,fraction=fraction,setIntersect=setIntersect);

		model.formula <- formula(bk$backfrm);
		nmodel <- modelFitting(model.formula,data,type);
		changes = as.integer(bk$Removed);
		weight <- 1.0;
		if ((changes>0) && (!inherits(nmodel, "try-error")))
		{
#			cat ("RMSE :",bk$BootModel$testRMSE,"\n")
			if (bk$testRMSE >= minbootRMSE) 
			{
				weight <- 1.0;
			}
			else
			{
				minbootRMSE <- bk$testRMSE;
				weight <- 0.0;	# suppress the betas beyond the minimum
			}
			if (bk$testRMSE > stopFratio*minbootRMSE) 
			{
				changes=0;	# stop with an increase from the minimum test error
			}
			if (changes>0)
			{
				loopsAux = loopsAux + 1;
				changed <- 1;
				changes2 <- as.character(as.list(attr(terms(model),"variables")))[which(!(as.character(as.list(attr(terms(model),"variables")))%in%as.character(as.list(attr(terms(nmodel),"variables")))))]
				if (length(changes2)>1)
				{
					changes2<-changes2[2]
				}
				if (adjsize>1)
				{
					if ((length(beforeFSCmodel$coefficients)>0)&&(length(model$coefficients)>0))
					{
						for (i in 1:length(beforeFSCmodel$coefficients))
						{
							notadded = TRUE;
							for (j in 1:length(nmodel$coefficients))
							{
								if (names(beforeFSCmodel$coefficients)[i] == names(nmodel$coefficients)[j])
								{
									beforeFSCmodel$coefficients[i] <- (weight*beforeFSCmodel$coefficients[i] + (1.0-weight)*nmodel$coefficients[j]); # it will weight based on the probability
									notadded=FALSE;
								}
							}
							if (notadded)
							{
								beforeFSCmodel$coefficients[i] <- weight*beforeFSCmodel$coefficients[i]; # it will average with zero
								wts[i] = wts[i]*weight;
							}
						}
					}
				}
				model <- nmodel;
			}
		}
		else
		{
			model <- nmodel;
			if (changes == 0) changes = 1;
			loopsAux = loopsAux + 1;      
		}
	}
	if (changes>=0)
	{
		modelReclas <- getVar.Res(model,data=data,Outcome=Outcome,type);
		NeRICV <- bootstrapValidation_Res(fraction,loops,model.formula,Outcome,data,type,plots=plots);
	}
	# else
	# {
		# NeRICV <-  NULL;
	# }
	if ((adjsize>1)&&(print==TRUE)) 
	{
		cat("Before BSC Mod:",beforeFSC.model.formula,"\n");
		if (!is.null(bk)) cat("Final  Formula:",bk$backfrm,"\n")
		cat("Adjust size:",adjsize,"\n");
		cat("Start RMSE:",startRMSE,"Min RMSE:",minbootRMSE,"final RMSE:",NeRICV$testRMSE,"\n")
	}


	result <- list(back.model=model,
	loops=loopsAux,
	reclas.info=modelReclas,
	bootCV=NeRICV,
	back.formula=formula(model.formula),
	lastRemoved=changes2,
	beforeFSC.model = beforeFSCmodel,
	beforeFSC.formula = formula(beforeFSC.model.formula));
	
	return (result);
}