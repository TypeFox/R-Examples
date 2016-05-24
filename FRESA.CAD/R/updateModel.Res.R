updateModel.Res <-
function(Outcome,covariates="1",pvalue=c(0.025,0.05),VarFrequencyTable,variableList,data,type=c("LM","LOGIT","COX"),testType=c("Binomial","Wilcox","tStudent"), lastTopVariable= 0,timeOutcome="Time",interaction=1,maxTrainModelSize=-1,bootLoops=1)
{
	type <- match.arg(type)

	acovariates <- covariates[1];
	if (length(covariates)>1)
	{
		for (i in 2:length(covariates))
		{	
			acovariates <- paste(acovariates," + ",covariates[i])
		}
	}
	covariates <- acovariates;

	varsize = ncol(data)-1;
	modsize=0;
	maxp <- max(pvalue); # set the pvalue to max so it will get the largest probable size of the non-corrected model
	loopst = 2;
	climpvalue = 0.4999 # 	 set to one third
	
	vnames <- as.vector(variableList[,1]);
	topvarID <- as.numeric(rownames(VarFrequencyTable));
	vnames_model <- vector();
	nsize <- nrow(data)
	if (maxTrainModelSize <= 0)
	{
# at least 5 samples per feature
		maxTrainModelSize = as.integer(nsize/5);
	}

	
	baseForm = Outcome;
#For Cox  models 
	if (type == "COX")
	{
	  baseForm = paste("Surv(",timeOutcome);
	  baseForm = paste(baseForm,paste(",",paste(Outcome,")")));
	}
	varlist <- vector();

	frm1 = paste(baseForm,paste(" ~ ",covariates));
	frm1 <- paste(frm1," + ");
	frm1 <- paste(frm1,vnames[topvarID[1]]);
	
	ftmp <- formula(frm1);
	varlist <- append(varlist,topvarID[1])
#	cat(frm1," start\n");
	bestmodel <- modelFitting(ftmp,data,type,TRUE)
	startIndex = 2;
	if ( inherits(bestmodel, "try-error"))
	{
		frm1 <- paste(frm1," + ",vnames[topvarID[2]]);
		varlist <- append(varlist,topvarID[2])
		ftmp <- formula(frm1);
		bestmodel <- modelFitting(ftmp,data,type,TRUE);
		startIndex = 3;
		topvarID[2]=0;
	}
	else
	{
		topvarID[1]=0;
	}
	
	bestResiduals <- residualForFRESA(bestmodel,data,Outcome);
	myTrainSample <- data;

	model_ziri <- vector();

	loops = 0;
	changes = 1;
	if (lastTopVariable < 1) lastTopVariable = length(VarFrequencyTable);
	if (lastTopVariable > length(VarFrequencyTable)) lastTopVariable = length(VarFrequencyTable);
	inserted = 1
	kins=1
	cpyformula <- frm1;
	termsinserted = 1;
#	cat("size ->",modsize,"\n")
	theBootLoops=bootLoops
	for (pval in 1:length(pvalue))
	{
		cthr_s = pvalue[pval];
		pthrO_s = cthr_s*cthr_s;
		loops = 0;
		changes = 1;
		if (pval>1) theBootLoops=1;
		while ((termsinserted < maxTrainModelSize)&&((changes>0) && (loops<loopst)))
		{
			tinserted = 0;
#			cthr_s = pvalue[pval]*(loops+1)/loopst;
#			cthr_s = pvalue[pval]/(loops+1);
			cthr_s = pvalue[pval];
			cthr = cthr_s;
			pthrO = min(pthrO_s,cthr^2);
#			if (changes > 0) cat ("Testing at :",cthr,"\n");

			changes = 0;
			samples <- nsize;
			
			ftmp <- formula(frm1);
#			cat("Update Formula 1: ",frm1,"\n")
			bestmodel <- modelFitting(ftmp,myTrainSample,type,TRUE)
			if ((loops == 0)&&(inherits(bestmodel, "try-error")))
			{
				frm1 <- paste(frm1," + ",vnames[topvarID[startIndex]]);
#				cat("Update Formula 1: ",frm1,"\n")
				varlist <- append(varlist,topvarID[startIndex]);
				VarFrequencyTable[startIndex]=0;
				ftmp <- formula(frm1);
				inserted = inserted + 1;
				tinserted = tinserted + 1;
				termsinserted = termsinserted + 1;
				bestmodel <- modelFitting(ftmp,myTrainSample,type,TRUE);
				startIndex = startIndex + 1;
			}
			
			while (inherits(bestmodel, "try-error"))
			{
				frm1 <- cpyformula;
				ftmp <- formula(frm1);
#				cat("Update Formula 2: ",frm1,"\n")

				bestmodel <- modelFitting(ftmp,myTrainSample,type,TRUE)
			}
			cpyformula <- frm1;

			bestResiduals <- residualForFRESA(bestmodel,myTrainSample,Outcome);

			for ( i in startIndex:lastTopVariable)
			{
				if ((VarFrequencyTable[i]>0) && (topvarID[i]>0) && (termsinserted < maxTrainModelSize))
				{
					frma <- paste(frm1," + ");
					frma <-paste(frma,vnames[topvarID[i]]);
					ftmp <- formula(frma);
					newmodel <- modelFitting(ftmp,myTrainSample,type,TRUE)
					if ( !inherits(newmodel, "try-error"))
					{
						if (bootLoops<4)
						{
							iprob <- .Call("improvedResidualsCpp",bestResiduals,residualForFRESA(newmodel,myTrainSample,Outcome),testType,0);
							piri <- iprob$p.value;
						}
						else
						{
#							cat("Update Formula 1: ",frma,"\n")
							iprob <- .Call("improvedResidualsCpp",bestResiduals,residualForFRESA(newmodel,myTrainSample,Outcome),testType,0);
							piri <- iprob$p.value;
							if (piri<cthr)
							{
#								cat("Update Formula 2: ",frma,"\n")
								bootmodel <- bootstrapValidation_Res(1.0000,theBootLoops,ftmp,Outcome,data,type,plots=FALSE)
								lastc <- ncol(bootmodel$tStudent.pvalues)
								
								switch(testType, 
									tStudent = 
									{ 
										ci <- as.vector(quantile(bootmodel$tStudent.pvlaues[,lastc], probs = c(climpvalue, 0.5,1.0 - climpvalue), na.rm = TRUE,names = FALSE, type = 7));
										ci2 <- median(bootmodel$test.tStudent.pvalues[,lastc]);
									},
									Wilcox = 
									{ 
										ci <- as.vector(quantile(bootmodel$wilcox.pvlaues[,lastc], probs = c(climpvalue, 0.5,1.0 - climpvalue), na.rm = TRUE,names = FALSE, type = 7));
										ci2 <- median(bootmodel$test.wilcox.pvalues[,lastc]);
									},
									Binomial =
									{ 
										ci <- as.vector(quantile(bootmodel$bin.pvlaues[,lastc], probs = c(climpvalue, 0.5,1.0 - climpvalue), na.rm = TRUE,names = FALSE, type = 7));
										ci2 <- median(bootmodel$test.bin.pvlaues[,lastc]);
									},
									Ftest =
									{ 
										ci <- as.vector(quantile(bootmodel$F.pvlaues[,lastc], probs = c(climpvalue, 0.5,1.0 - climpvalue), na.rm = TRUE,names = FALSE, type = 7));
										ci2 <- median(bootmodel$test.F.pvlaues[,lastc]);
									},
								)
								piri <- max(ci[3],ci2,iprob$p.value)
#								piri <- max(ci[3],ci2,iprob$p.value); # the maximum of the train, test median or sample p-value
#								cat(vnames[topvarID[i]],"train: ",ci," Test: ",ci2," All: ",iprob$p.value," Max: ",piri,"\n")
							}
						}


						if (is.numeric(piri) && !is.na(piri) && (piri<cthr))
						{
							bestResiduals <- residualForFRESA(newmodel,myTrainSample,Outcome);
							frm1 <- paste(frm1," + ",vnames[topvarID[i]]);
							vnames_model <- append(vnames_model,vnames[topvarID[i]]);
							varlist <- append(varlist,topvarID[i]);
							model_ziri <- append(model_ziri,abs(qnorm(piri)));
							changes = changes + 1;
							inserted = inserted + 1;
							tinserted = tinserted + 1;
							termsinserted = termsinserted + 1;
							kins=1
							VarFrequencyTable[i]=0;
						}	
						if (interaction == 2)
						{
							for (nlist in 1:inserted)
							{
								if (termsinserted < maxTrainModelSize)
								{
									if (kins==1)
									{
										pthrOl=cthr;
										frma <- paste(frm1," + I(",vnames[varlist[nlist]],"*",vnames[topvarID[i]],")")
									}
									else
									{
										frma <- paste(frm1," + ",vnames[topvarID[i]]," + I(",vnames[varlist[nlist]],"*",vnames[topvarID[i]],")")
										pthrOl=pthrO;
									}
									ftmp <- formula(frma);
									newmodel <- modelFitting(ftmp,myTrainSample,type,TRUE)
									if ( !inherits(newmodel, "try-error"))
									{
										if (bootLoops<4)
										{
											iprob <- .Call("improvedResidualsCpp",bestResiduals,residualForFRESA(newmodel,myTrainSample,Outcome),testType,0);
											piri <- iprob$p.value;
										}
										else
										{
											iprob <- .Call("improvedResidualsCpp",bestResiduals,residualForFRESA(newmodel,myTrainSample,Outcome),testType,0);
											piri <- iprob$p.value;
											if (piri<pthrOl)
											{
												bootmodel <- bootstrapValidation_Res(1.0000,theBootLoops,ftmp,Outcome,data,type,plots=FALSE)
												lastc <- ncol(bootmodel$tStudent.pvalues)
												switch(testType, 
													tStudent = 
													{ 
														ci <- as.vector(quantile(bootmodel$tStudent.pvlaues[,lastc], probs = c(climpvalue, 0.5,1.0 - climpvalue), na.rm = TRUE,names = FALSE, type = 7));
														ci2 <- median(bootmodel$test.tStudent.pvalues[,lastc]);
													},
													Wilcox = 
													{ 
														ci <- as.vector(quantile(bootmodel$wilcox.pvlaues[,lastc], probs = c(climpvalue, 0.5,1.0 - climpvalue), na.rm = TRUE,names = FALSE, type = 7));
														ci2 <- median(bootmodel$test.wilcox.pvalues[,lastc]);
													},
													Binomial =
													{ 
														ci <- as.vector(quantile(bootmodel$bin.pvlaues[,lastc], probs = c(climpvalue, 0.5,1.0 - climpvalue), na.rm = TRUE,names = FALSE, type = 7));
														ci2 <- median(bootmodel$test.bin.pvlaues[,lastc]);
													},
													Ftest =
													{ 
														ci <- as.vector(quantile(bootmodel$F.pvlaues[,lastc], probs = c(climpvalue, 0.5,1.0 - climpvalue), na.rm = TRUE,names = FALSE, type = 7));
														ci2 <- median(bootmodel$test.F.pvlaues[,lastc]);
													},
												)
												piri <- max(ci[3],ci2,iprob$p.value);
											}
										}
										if (is.numeric(piri) && !is.na(piri) && (piri<pthrOl))
										{
											bestResiduals <- residualForFRESA(newmodel,myTrainSample,Outcome);
											frm1 <- frma;
											vnames_model <- append(vnames_model,vnames[topvarID[i]]);
											model_ziri <- append(model_ziri,abs(qnorm(piri)));
											if (kins == 0)
											{
												varlist <- append(varlist,topvarID[i]);
												inserted = inserted + 1;
												tinserted = tinserted + 1;
												termsinserted = termsinserted + 1;
											}
											tinserted = tinserted + 1;
											termsinserted = termsinserted + 1;
											kins =1
											VarFrequencyTable[i]=0;
											changes = changes + 1;
										}
									}
								}
							}
						}
					}
					kins=0
				}
			}
			loops = loops+1;
		}
#		cat(frm1,"\n");
	}

	ftmp <- formula(frm1);
	bestmodel <- modelFitting(ftmp,data,type)
#	cat("Update Formula:",frm1,"\n");
#	print(summary(bestmodel));
	
  	result <- list(final.model=bestmodel,
	var.names=vnames_model,
	formula=ftmp,
	z.NeRI=model_ziri,
	loops=loops);
  
	return (result);
}
