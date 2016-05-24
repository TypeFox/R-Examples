updateModel.Bin <-
function(Outcome,covariates="1",pvalue=c(0.025,0.05),VarFrequencyTable,variableList,data,type=c("LM","LOGIT","COX"), lastTopVariable= 0,timeOutcome="Time",selectionType=c("zIDI","zNRI"),numberOfModels=3,interaction=1,maxTrainModelSize=0,bootLoops=1)
{
	type <- match.arg(type)
  	seltype <- match.arg(selectionType)

	acovariates <- covariates[1];
	if (length(covariates)>1)
	{
		for (i in 2:length(covariates))
		{	
			acovariates <- paste(acovariates," + ",covariates[i])
		}
	}
	covariates <- acovariates;

	
	vnames <- as.vector(variableList[,1]);
	topvarID <- as.numeric(rownames(VarFrequencyTable));
	
	casesample = subset(data,get(Outcome)  == 1);
	controlsample = subset(data,get(Outcome) == 0);

	sizecases = nrow(casesample);
	sizecontrol = nrow(controlsample);
	minsize = min(sizecases,sizecontrol);
	varsize = ncol(data)-1;

	if (numberOfModels <= 0) 
	{
		numberOfModels = 1;
	}
	if (maxTrainModelSize == 0)
	{
		maxTrainModelSize = as.integer(minsize/5);
	}

	
	baseForm = Outcome;
					
	#For Cox  models 
	if (type == "COX")
	{
		baseForm = paste("Surv(",timeOutcome);
		baseForm = paste(baseForm,paste(",",paste(Outcome,")")));
	}
		
	vnames_model <- vector();
	model_zmin <- vector();
	topvarCpy <- topvarID

	
	if (lastTopVariable < 1) lastTopVariable = length(VarFrequencyTable);
	if (lastTopVariable > length(VarFrequencyTable)) lastTopVariable = length(VarFrequencyTable);
	firstVar = 1
	formulaList <- vector();
	varlist <- vector();
	inserted = 0
	loops = 0
	bestmodel <- NULL;
	ftmp <- NULL;
	frm1 <- NULL;
	modsize = 0;
	maxp <- max(pvalue);	# set the pvalue to maxp so it will get the most probable size of the non-corrected model 
	
	mysample <- data;
	climpvalue = 0.4999;	# 
	loopst = 2;
	maxp <- min(2.0*maxp,0.20);
	theBootLoops=bootLoops;
	if (lastTopVariable>1)
	{


		for (nMod in 1:numberOfModels)
		{
			loops = 0;
			inserted = 0;
			termsinserted = 0;
			if (firstVar>0)
			{
				varlist <- vector();
				for (pidx in 1:length(pvalue))
				{
					cthr = cthr_s = abs(qnorm(pvalue[pidx]));	
					zthrO = zthrO_s = abs(qnorm((pvalue[pidx]^2)));
					changes = 1;
					if (pidx>1) theBootLoops=1;
					while ((termsinserted < maxTrainModelSize)&&((changes>0) && (loops<loopst)))
					{
#						cthr_s = abs(qnorm(pvalue[pidx]*(loops+1)))/loopst;
#						cthr_s = abs(qnorm(pvalue[pidx]/(loops+1)));
						cthr_s = abs(qnorm(pvalue[pidx]));
						zthrO_s = abs(qnorm((pvalue[pidx]^2)));
						cthr = cthr_s;
						zthrO = zthrO_s;
						tinserted = 0;
#						if (changes > 0) cat ("Testing at :",cthr,"\n");
						changes = 0;

						if ((loops == 0)&&(topvarID[firstVar]>0))
						{
							frm1 = paste(baseForm,paste(" ~ ",covariates));
							frm1 <- paste(frm1," + ");
							frm1 <- paste(frm1,vnames[topvarID[firstVar]]);
							varlist <- append(varlist,topvarID[firstVar]);
							vnames_model <- append(vnames_model,vnames[topvarID[firstVar]]);
							model_zmin <- append(model_zmin,NA);
							topvarID[firstVar] = 0;
							inserted = 1;
							tinserted = 1;
							termsinserted = termsinserted + 1;
							ftmp <- formula(frm1);
							bestmodel <- modelFitting(ftmp,mysample,type,TRUE)
						}
						else
						{									
							ftmp <- formula(frm1);
							bestmodel <- modelFitting(ftmp,mysample,type,TRUE)							
						}
						sizetrain <- nrow(mysample);
						if ( !inherits(bestmodel, "try-error"))
						{
							bestpredict <- predictForFresa(bestmodel,mysample,'prob');

							firstVar = 1;
							for ( i in 2:lastTopVariable)
							{
								if ((VarFrequencyTable[i]>0) && (topvarID[i]>0) && (termsinserted < maxTrainModelSize))
								{

									kinserted = 0
									frma <- paste(frm1," + ");
									frma <-paste(frma,vnames[topvarID[i]]);
									ftmp <- formula(frma);
									newmodel <- modelFitting(ftmp,mysample,type,TRUE)
									if ( !inherits(newmodel, "try-error"))
									{
										iprob_t <- .Call("improveProbCpp",bestpredict,predictForFresa(newmodel,mysample,'prob'),mysample[,Outcome],0);
										if (seltype=="zIDI") 
										{
											zmin <- iprob_t$z.idi;
										}
										else
										{
											zmin <- iprob_t$z.nri;
										}
										if (is.numeric(zmin) && !is.na(zmin))
										{
											if ((theBootLoops>4)&&(zmin>cthr))
											{
												idiCV <- bootstrapValidation_Bin(1.0000,theBootLoops,ftmp,Outcome,mysample,type,plots=FALSE)
												lastc <- ncol(idiCV$z.IDIs)
												if (seltype=="zIDI")
												{
													ci <- as.vector(quantile(idiCV$z.IDIs[,lastc], probs = c(climpvalue, 0.5, 0.75), na.rm = TRUE,names = FALSE, type = 7));
													ci2 <- median(idiCV$test.z.IDIs[,lastc]);
												}
												else
												{
													ci <- as.vector(quantile(idiCV$z.NRIs[,lastc], probs = c(climpvalue, 0.5, 0.75), na.rm = TRUE,names = FALSE, type = 7));
													ci2 <- median(idiCV$test.z.NRIs[,lastc]);
												}
												zmin = min(ci[1],ci2,zmin);
											}

											if (is.numeric(zmin) && !is.na(zmin)&&(zmin>cthr))
											{
												bestpredict <-predictForFresa(newmodel,mysample,'prob');

												frm1 <- frma;
												vnames_model <- append(vnames_model,vnames[topvarID[i]]);
												model_zmin <- append(model_zmin,zmin);
												varlist <- append(varlist,topvarCpy[i]);
												changes = changes + 1;
												inserted = inserted + 1;
												tinserted = tinserted + 1;
												termsinserted = termsinserted + 1;
												topvarID[i] = 0
											}
											if (is.numeric(zmin) && !is.na(zmin)&&(zmin<=cthr))
											{
												if (firstVar == 1) firstVar = i;
											}
										}
										if (interaction == 2)
										{
											chkin <- (topvarID[i] > 0);
											for (nlist in 1:inserted)
											{
												if (termsinserted < maxTrainModelSize)
												{
													if (topvarID[i] == 0)
													{
														frma <- paste(frm1," + I(",vnames[varlist[nlist]],"*",vnames[topvarCpy[i]],")")
														zthrOl = cthr;
													}
													else
													{
														frma <- paste(frm1," + ",vnames[topvarCpy[i]]," + I(",vnames[varlist[nlist]],"*",vnames[topvarCpy[i]],")")
														zthrOl = zthrO;
													}
													ftmp <- formula(frma);
													newmodel <- modelFitting(ftmp,mysample,type,TRUE)
													if ( !inherits(newmodel, "try-error"))
													{
														iprob_t <- .Call("improveProbCpp",bestpredict,predictForFresa(newmodel,mysample,'prob'),mysample[,Outcome],0);
														if (seltype=="zIDI") 
														{
															zmin <- iprob_t$z.idi;
														}
														else
														{
															zmin <- iprob_t$z.nri;
														}
														if (is.numeric(zmin)&&!is.na(zmin))
														{
															if ((bootLoops>4)&&(zmin>zthrOl))
															{
																idiCV <- bootstrapValidation_Bin(1.0000,theBootLoops,ftmp,Outcome,mysample,type,plots=FALSE)
																lastc <- ncol(idiCV$z.IDIs)
																if (seltype=="zIDI")
																{
																	ci <- as.vector(quantile(idiCV$z.IDIs[,lastc], probs = c(climpvalue, 0.5, 0.75), na.rm = TRUE,names = FALSE, type = 7));
																	ci2 <- median(idiCV$test.z.IDIs[,lastc]);
																}
																else
																{
																	ci <- as.vector(quantile(idiCV$z.NRIs[,lastc], probs = c(climpvalue, 0.5, 0.75), na.rm = TRUE,names = FALSE, type = 7));
																	ci2 <- median(idiCV$test.z.NRIs[,lastc]);
																}
																zmin = min(ci[1],ci2,zmin);
															}

															if (is.numeric(zmin)&&!is.na(zmin)&&(zmin>zthrOl))
															{
																bestpredict <- predictForFresa(newmodel,mysample,'prob');
																
																frm1 <- frma;
																vnames_model <- append(vnames_model,vnames[topvarCpy[i]]);
																model_zmin <- append(model_zmin,zmin);
																changes = changes + 1;
																kinserted = kinserted + 1;
																termsinserted = termsinserted + 1;
																tinserted = tinserted + 1;
																if (topvarID[i] != 0)
																{
																	termsinserted = termsinserted + 1;
																	tinserted = tinserted + 1;
																}
																topvarID[i] = 0
															}
														}
													}
												}
											}
											if ((kinserted > 0) && chkin ) 
											{
												varlist <- append(varlist,topvarCpy[i]);
												inserted = inserted + 1;
												if (firstVar == 1) firstVar = i;
											}
										}
									}
								}
							}
							loops = loops+1;
						}
					}
#					cat (frm1,"\n")
				}
				formulaList <- append(formulaList,frm1)
				ftmp <- formula(frm1);
			}
		}
	}
	if (length(formulaList)==0)
	{
		frm1 = paste(baseForm," ~ ",covariates," + ",vnames[topvarID[1]]);
#		cat ("Update Formula: ",frm1,"\n")
		formulaList <- append(formulaList,frm1)
		ftmp <- formula(formulaList[1]);
		bestmodel <- modelFitting(ftmp,data,type)
	}
	else
	{
#		cat("Top Formula: \n")
		ftmp <- formula(formulaList[1]);
#		cat ("Update Formula: ",formulaList[1],"\n")
		bestmodel <- modelFitting(ftmp,data,type)
	}
#	print(summary(bestmodel));

#	cat ("Update Formula:",formulaList[1],"\n")
	
  	result <- list(final.model=bestmodel,
	var.names=vnames_model,
	formula=ftmp,
	z.selectionType=model_zmin,
	loops=loops,
	formula.list=formulaList);
  
	return (result);
}
