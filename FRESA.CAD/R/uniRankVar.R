uniRankVar <-
function(variableList,formula,Outcome,data,categorizationType=c("Raw","Categorical","ZCategorical","RawZCategorical","RawTail","RawZTail"),type=c("LOGIT","LM","COX"),rankingTest=c("zIDI","zNRI","IDI","NRI","NeRI","Ztest","AUC","CStat","Kendall"),cateGroups=c(0.1,0.9),raw.dataFrame=NULL,description=".",uniType=c("Binary","Regression"),FullAnalysis=TRUE) 
{


	if (is.null(raw.dataFrame))  raw.dataFrame <- data;
	type <- match.arg(type);
	uniType <- match.arg(uniType);
	categorizationType <- match.arg(categorizationType);
	rankingTest <- match.arg(rankingTest);
	colnamesList <- as.vector(variableList[,1]);


	if (description == ".")
	{
		descripList <- colnamesList;
	}
	else
	{
		descripList <- as.vector(variableList[,description]);
	}
	
	if (uniType=="Binary")
	{
		caserawsample <- subset(raw.dataFrame,get(Outcome)  == 1);
		controlrawsample <- subset(raw.dataFrame,get(Outcome) == 0);

		caseZsample <- subset(data,get(Outcome)  == 1);
		controlZsample <- subset(data,get(Outcome) == 0);

		sizecaseZsample <- nrow(caseZsample);
		sizecontrolZsample <- nrow(controlZsample);
	}
	else
	{
		caserawsample <- NULL;
		controlrawsample <- NULL;

		caseZsample <- NULL;
		controlZsample <- NULL;

		sizecaseZsample <- NULL;
		sizecontrolZsample <- NULL;
	}

	size = length(colnamesList);
	Name <- vector();

	parent <- vector();
	descrip <- vector();


	IDI <- vector();
	NRI <- vector();
	zIDI <- vector();
	zNRI <- vector();
	ROCAUC <- vector();
	ZGLM <- vector();

	NeRI <- vector();
	BinRes.p <- vector();
	WilcoxRes.p <- vector();
	TstudentRes.p <- vector();
	FRes.p <- vector();


	cohortMean <- vector();
	cohortStd <- vector();
	cohortKSD <- vector();
	cohortKSP <- vector();
	cohortZKSP <- vector();
	cohortZKSD <- vector();

	
	caseMean <- vector();
	caseStd <- vector();
	caseKSD <- vector();
	caseKSP <- vector();
	caseZKSP <- vector();
	caseZKSD <- vector();
	caseN_Z_Low_Tail <- vector();
	caseN_Z_Hi_Tail <- vector();

	controlMean <- vector();
	controlStd <- vector();
	controlKSD <- vector();
	controlKSP <- vector();
	controlZKSP <- vector();
	controlZKSD <- vector();
	controlN_Z_Low_Tail <- vector();
	controlN_Z_Hi_Tail <- vector();
	kendall.r <- vector();
	spearman.r <- vector();
	pearson.r <- vector();
	kendall.p <- vector();
	cStatCorr <- vector();
	
	t.Rawvalue <- vector();
	t.Zvalue <- vector();
	wilcox.Zvalue <- vector();
	
	frm1 <- formula;
	ftmp <- formula(frm1);
	bmodel <- modelFitting(ftmp,data,type)
	baseResiduals <- residualForFRESA(bmodel,data,Outcome)+rnorm(nrow(data),0,1e-10);
	basepredict <- predictForFresa(bmodel,data, 'prob');
	if (type=="COX")
	{
		termslist <- attr(terms(ftmp),"term.labels");
		if (length(termslist)==0)
		{
			ftmp <- formula(paste(frm1,"+dummy"));
			dummy <-  rnorm(nrow(data));
			coxframe <- cbind(data,dummy);
			bmodel <- modelFitting(ftmp,coxframe,type)
			basepredict <- predictForFresa(bmodel,coxframe, 'prob');
			baseResiduals <- residualForFRESA(bmodel,coxframe,Outcome);
		}
	}
	
#	if (!FullAnalysis)
#	{
#		rankingTest="Ztest";
#	}
	
	for (j in 1:size)
	{
#		cat (colnamesList[j],"\n")
		frm1 = formula;
		categories = 1;
		catlist <- vector();
		caseCount1 = 0;
		caseCount2 = 0;
		controlCount1 = 0;
		controlCount2 = 0;
		stddf=0;
		kendcor <- NA;
		pearcor <- NA;
		speacor <- NA;
		cstat <- NA;
		kstZdf <- NA;
		kstdf <- NA;
		meCa <- NA;
		stdCa <- NA;
		kstCa <- NA;
		kstZCa <- NA;
		meCo <- NA;
		stdCo <- NA;
		kstCo <- NA;
		kstZCo <- NA;
		rtt <- NA;
		ztt <- NA;
		medf <- NA;
		if ((uniType=="Binary")&& FullAnalysis)
		{
			wtt <- -qnorm(wilcox.test(controlrawsample[,colnamesList[j]],caserawsample[,colnamesList[j]],na.action=na.exclude)$p.value,0,1);
		}
		else
		{
			wtt <- NA;
		}
#		cat (colnamesList[j],"Wilcox: ",wtt,"\n")
		if (FullAnalysis)
		{
			stddf <- sd(raw.dataFrame[,colnamesList[j]],na.rm = TRUE);
			if (stddf>0)
			{
				kendcor <- try(cor.test(data[,colnamesList[j]],data[,Outcome],method="kendall",na.action=na.omit));
				pearcor <- try(cor.test(data[,colnamesList[j]],data[,Outcome],method="pearson",na.action=na.omit));
				speacor <- try(cor.test(data[,colnamesList[j]],data[,Outcome],method="spearman",na.action=na.omit));
				if (inherits(kendcor, "try-error")) kendcor$estimate = 0;
				if (is.na(kendcor$estimate)) kendcor$estimate = 0;
				if (kendcor$estimate > 0)
				{
					cstat <- rcorr.cens(data[,colnamesList[j]],data[,Outcome], outx=FALSE)
				}
				else
				{
					cstat <- rcorr.cens(data[,colnamesList[j]],-data[,Outcome], outx=FALSE)
				}
			}
		}
#		cat (colnamesList[j],"cstat: ",cstat,"\n")
		if (length(table(data[,colnamesList[j]]))>4)
		{

			if (FullAnalysis)
			{
				medf <- mean(data[,colnamesList[j]],na.rm = TRUE);
				stddf <- sd(data[,colnamesList[j]],na.rm = TRUE);
				kstZdf <- ks.test(data[,colnamesList[j]],"pnorm",medf,stddf);
				medf <- mean(raw.dataFrame[,colnamesList[j]],na.rm = TRUE);
				stddf <- sd(raw.dataFrame[,colnamesList[j]],na.rm = TRUE);
				kstdf <- ks.test(raw.dataFrame[,colnamesList[j]],"pnorm",medf,stddf);
			}
			if (uniType=="Binary")
			{
				if (FullAnalysis)
				{
					meCa <- mean(caseZsample[,colnamesList[j]],na.rm = TRUE);
					stdCa <- sd(caseZsample[,colnamesList[j]],na.rm = TRUE);
					kstZCa <- ks.test(caseZsample[,colnamesList[j]],"pnorm",meCa,stdCa);
					
					meCo <- mean(controlZsample[,colnamesList[j]],na.rm = TRUE);
					stdCo <- sd(controlZsample[,colnamesList[j]],na.rm = TRUE);
					kstZCo <- ks.test(controlZsample[,colnamesList[j]],"pnorm",meCo,stdCo);
					
					meCa <- mean(caserawsample[,colnamesList[j]],na.rm = TRUE);
					stdCa <- sd(caserawsample[,colnamesList[j]],na.rm = TRUE);
					kstCa <- ks.test(caserawsample[,colnamesList[j]],"pnorm",meCa,stdCa);
					
					meCo <- mean(controlrawsample[,colnamesList[j]],na.rm = TRUE);
					stdCo <- sd(controlrawsample[,colnamesList[j]],na.rm = TRUE);
					kstCo <- ks.test(controlrawsample[,colnamesList[j]],"pnorm",meCo,stdCo);
					
					rtt <- try(t.test(controlrawsample[,colnamesList[j]],caserawsample[,colnamesList[j]],na.action=na.exclude));
					ztt <- try(t.test(controlZsample[,colnamesList[j]],caseZsample[,colnamesList[j]],na.action=na.exclude));
				}

				if (!is.na(cateGroups[1]))
				{
					if (categorizationType != "Raw")
					{
						zthr = sprintf("[,'%s'] < %5.3f )",colnamesList[j],qnorm(cateGroups[1]));

						caseCount1 <- eval(parse(text = paste("sum(caseZsample",zthr)));
						controlCount1 <- eval(parse(text = paste("sum(controlZsample",zthr)));

						categories=length(cateGroups);
						if (!is.na(cateGroups[categories]))
						{
							zthr = sprintf("[,'%s'] < %5.3f )",colnamesList[j],qnorm(cateGroups[categories]));
							caseCount2 <- eval(parse(text = paste("sum(caseZsample",zthr)));
							controlCount2 <- eval(parse(text = paste("sum(controlZsample",zthr)));

						}
						else
						{
							zthr = sprintf("[,'%s'] < %5.3f )",colnamesList[j],1-qnorm(cateGroups[1]));
							caseCount2 <- eval(parse(text = paste("sum(caseZsample",zthr)));
							controlCount2 <- eval(parse(text = paste("sum(controlZsample",zthr)));
						}
					}
				}
			}


			if ((min(data[,colnamesList[j]])<0) && (length(table(data[,colnamesList[j]]))>4))
			{
				switch(categorizationType,
					Raw =
					{
						categories=1;
						catlist <- append(catlist,colnamesList[j]);				
					},
					Categorical =
					{
						categories=length(cateGroups);

						zthr = sprintf("%5.3f",qnorm(cateGroups[1]));


						for (n in 1:categories)
						{
							if (n==1)
							{
								zthr = sprintf("%5.3f",qnorm(cateGroups[n]));
								catvar = paste("I(",colnamesList[j]);
								catvar = paste(catvar," < ");
								catvar = paste(catvar,zthr);
								catvar = paste(catvar,")");
								catlist <- append(catlist,catvar);
							}
							else
							{
								zthr = sprintf("%5.3f",qnorm(cateGroups[n-1]));
								zthr2 = sprintf("%5.3f",qnorm(cateGroups[n]));
								catvar = paste("I((",colnamesList[j]);
								catvar = paste(catvar," >= ");
								catvar = paste(catvar,zthr);
								catvar = paste(catvar,") & (");
								catvar = paste(catvar,colnamesList[j]);
								catvar = paste(catvar," < ");
								catvar = paste(catvar,zthr2);
								catvar = paste(catvar,"))");
								catlist <- append(catlist,catvar);
							}
						}
						zthr = sprintf("%5.3f",qnorm(cateGroups[categories]));
						catvar = paste("I(",colnamesList[j]);
						catvar = paste(catvar," >= ");
						catvar = paste(catvar,zthr);
						catvar = paste(catvar,")");
						catlist <- append(catlist,catvar);
						categories = categories + 1;
					},
					ZCategorical =
					{
						categories=length(cateGroups);
						zthr = sprintf("%5.3f",qnorm(cateGroups[1]));

						for (n in 1:categories)
						{

							if (n==1)
							{
								zthr = sprintf("%5.3f",qnorm(cateGroups[n]));
								catvar = paste("I(",colnamesList[j]);
								catvar = paste(catvar,"* (");
								catvar = paste(catvar,colnamesList[j]);
								catvar = paste(catvar," < ");
								catvar = paste(catvar,zthr);
								catvar = paste(catvar,"))");
								catlist <- append(catlist,catvar);
							}
							else
							{
								zthr = sprintf("%5.3f",qnorm(cateGroups[n-1]));
								zthr2 = sprintf("%5.3f",qnorm(cateGroups[n]));
								catvar = paste("I(",colnamesList[j]);
								catvar = paste(catvar,"* ((");
								catvar = paste(catvar,colnamesList[j]);
								catvar = paste(catvar," >= ");
								catvar = paste(catvar,zthr);
								catvar = paste(catvar,") & (");
								catvar = paste(catvar,colnamesList[j]);
								catvar = paste(catvar," < ");
								catvar = paste(catvar,zthr2);
								catvar = paste(catvar,")))");
								catlist <- append(catlist,catvar);
							}
						}
						zthr = sprintf("%5.3f",qnorm(cateGroups[categories]));
						catvar = paste("I(",colnamesList[j]);
						catvar = paste(catvar,"* (");
						catvar = paste(catvar,colnamesList[j]);
						catvar = paste(catvar," >= ");
						catvar = paste(catvar,zthr);
						catvar = paste(catvar,"))");
						catlist <- append(catlist,catvar);
						categories = categories+1;
					},
					RawZCategorical =
					{
						categories=length(cateGroups);

						zthr = sprintf("%5.3f",qnorm(cateGroups[1]));

						catlist <- append(catlist,colnamesList[j]);				
						for (n in 1:categories)
						{

							if (n==1)
							{
								zthr = sprintf("%5.3f",qnorm(cateGroups[n]));
								catvar = paste("I(",colnamesList[j]);
								catvar = paste(catvar,"* (");
								catvar = paste(catvar,colnamesList[j]);
								catvar = paste(catvar," < ");
								catvar = paste(catvar,zthr);
								catvar = paste(catvar,"))");
								catlist <- append(catlist,catvar);
							}
							else
							{
								zthr = sprintf("%5.3f",qnorm(cateGroups[n-1]));
								zthr2 = sprintf("%5.3f",qnorm(cateGroups[n]));
								catvar = paste("I(",colnamesList[j]);
								catvar = paste(catvar,"* ((");
								catvar = paste(catvar,colnamesList[j]);
								catvar = paste(catvar," >= ");
								catvar = paste(catvar,zthr);
								catvar = paste(catvar,") & (");
								catvar = paste(catvar,colnamesList[j]);
								catvar = paste(catvar," < ");
								catvar = paste(catvar,zthr2);
								catvar = paste(catvar,")))");
								catlist <- append(catlist,catvar);
							}
						}
						zthr = sprintf("%5.3f",qnorm(cateGroups[categories]));
						catvar = paste("I(",colnamesList[j]);
						catvar = paste(catvar,"* (");
						catvar = paste(catvar,colnamesList[j]);
						catvar = paste(catvar," >= ");
						catvar = paste(catvar,zthr);
						catvar = paste(catvar,"))");
						catlist <- append(catlist,catvar);
						categories = categories+2;
					},			
					RawZTail =
					{
						categories = 1;
						catlist <- append(catlist,colnamesList[j]);				

						if (!is.null(sizecaseZsample))
						{
							zthr = sprintf("%5.3f",qnorm(cateGroups[1]));
							f1= caseCount1/sizecaseZsample;
							f2= controlCount1/sizecontrolZsample;
							if ((f1>f2)&&(f1>0.1))
							{
								catvar = paste("I(",colnamesList[j]);
								catvar = paste(catvar,"* (");
								catvar = paste(catvar,colnamesList[j]);
								catvar = paste(catvar," < ");
								catvar = paste(catvar,zthr);
								catvar = paste(catvar,"))");
								catlist <- append(catlist,catvar);
								categories = categories+1;
							}				


							zthr = sprintf("%5.3f",qnorm(1.0-cateGroups[1]));
							f1= caseCount2/sizecaseZsample;
							f2= controlCount2/sizecontrolZsample;
							if ((f1>f2)&&(f1>0.1))
							{
								catvar = paste("I(",colnamesList[j]);
								catvar = paste(catvar,"* (");
								catvar = paste(catvar,colnamesList[j]);
								catvar = paste(catvar," > ");
								catvar = paste(catvar,zthr);
								catvar = paste(catvar,"))");
								catlist <- append(catlist,catvar);
								categories = categories+1;
							}
						}
					},			
					RawTail =
					{
						categories = 1;
						catlist <- append(catlist,colnamesList[j]);				
						if (!is.null(sizecaseZsample))
						{
							zthr = sprintf("%5.3f",qnorm(cateGroups[1]));
							f1 = caseCount1/sizecaseZsample;
							f2 = controlCount1/sizecontrolZsample; 
							if ((f1>f2)&&(f1>0.1))   # will add only if fraction is greater
							{
								catvar = paste("I(",colnamesList[j]);
								catvar = paste(catvar," < ");
								catvar = paste(catvar,zthr);
								catvar = paste(catvar,")");
								catlist <- append(catlist,catvar);
								categories = categories+1;
							}				


							zthr = sprintf("%5.3f",qnorm(1.0-cateGroups[1]));
							f1= caseCount2/sizecaseZsample;
							f2= controlCount2/sizecontrolZsample;
							if ((f1>f2)&&(f1>0.1)) # will add only if fraction is greater
							{
								catvar = paste("I(",colnamesList[j]);
								catvar = paste(catvar," > ");
								catvar = paste(catvar,zthr);
								catvar = paste(catvar,")");
								catlist <- append(catlist,catvar);
								categories = categories+1;
							}
						}
					},			
					{
						categories=1;
						catlist <- append(catlist,colnamesList[j]);
					}
				)
			}
			else
			{
				categories=1;
				catlist <- append(catlist,colnamesList[j]);
			}
		}
		else
		{
			categories=1;
			catlist <- append(catlist,colnamesList[j]);	
			medf = table(data[,colnamesList[j]])[1];			
			stddf = table(data[,colnamesList[j]])[2];			
			if (uniType=="Binary")
			{
				meCa <- table(caseZsample[,colnamesList[j]])[1];
				stdCa <- table(caseZsample[,colnamesList[j]])[2];
				
				meCo <- table(controlZsample[,colnamesList[j]])[1];
				stdCo <- table(controlZsample[,colnamesList[j]])[2];
			}
		}
		for (n in 1:categories)
		{
			termName <- str_replace_all(catlist[n]," ","");
			termName <- str_replace_all(termName,"<"," < ");
			termName <- str_replace_all(termName,">"," > ");
			termName <- str_replace_all(termName,"&"," & ");
			termName <- str_replace_all(termName,"=","= ");
			termName <- str_replace_all(termName,fixed("> ="),">=");
			termName <- str_replace_all(termName,fixed("*")," * ");
			frmg <- paste( formula,paste(" + ",termName));
#			cat(frmg,"\n")
			ftmg <- formula(frmg);
			if (type=="COX") 
			{
				zcol=4;
			}
			else
			{
				zcol=3;
			}
			lmodel <- modelFitting(ftmg,data,type)
#			print(summary(lmodel))
			if (!inherits(lmodel, "try-error"))
			{
				modcoef <- summary(lmodel)$coefficients;
				sizecoef <- length(lmodel$coef);
			}
			else
			{
				modcoef <- NULL;
				sizecoef <- NULL;
			}
				Name <- append(Name,termName);
				parent <- append(parent,colnamesList[j])
				descrip <- append(descrip,descripList[j])
				if (uniType=="Binary")
				{
					caseMean <- append(caseMean,meCa); 
					caseStd <- append(caseStd,stdCa);
					if (!is.na(kstCa[[1]])) 
					{
						caseKSD <- append(caseKSD,kstCa$statistic);
						caseKSP <- append(caseKSP,kstCa$p.value);
						caseZKSP <- append(caseZKSP,kstZCa$p.value);
						caseZKSD <- append(caseZKSD,kstZCa$statistic);
					}
					else
					{
						caseKSD <- append(caseKSD,NA);
						caseKSP <- append(caseKSP,NA);
						caseZKSP <- append(caseZKSP,NA);
						caseZKSD <- append(caseZKSD,NA);
					}
					controlMean <- append(controlMean,meCo); 
					controlStd <- append(controlStd,stdCo);
					if (!is.na(kstCo[[1]]))
					{
						controlKSD <- append(controlKSD,kstCo$statistic);
						controlKSP <- append(controlKSP,kstCo$p.value);
						controlZKSP <- append(controlZKSP,kstZCo$p.value);
						controlZKSD <- append(controlZKSD,kstZCo$statistic);
					}
					else
					{
						controlKSD <- append(controlKSD,NA);
						controlKSP <- append(controlKSP,NA);
						controlZKSP <- append(controlZKSP,NA);
						controlZKSD <- append(controlZKSD,NA);
					}
					if (!is.na(rtt[[1]]))
					{
						if ( !inherits(rtt, "try-error"))
						{
							t.Rawvalue <- append(t.Rawvalue,rtt$statistic);
						}
						else
						{
							t.Rawvalue <- append(t.Rawvalue,NA);
						}
						if ( !inherits(ztt, "try-error"))
						{
							t.Zvalue <- append(t.Zvalue,ztt$statistic);
						}
						else
						{
							t.Zvalue <- append(t.Zvalue,NA);
						}
					}
					else
					{
						t.Rawvalue <- append(t.Rawvalue,NA);
						t.Zvalue <- append(t.Zvalue,NA);
					}
					if (!is.na(wtt[[1]]))
					{
						wilcox.Zvalue <- append(wilcox.Zvalue,wtt);
					}
					else
					{
						wilcox.Zvalue <- append(wilcox.Zvalue,NA);
					}
				}
			
			cohortMean <- append(cohortMean,medf); 
			cohortStd <- append(cohortStd,stddf);
			if (!is.na(kstdf[[1]]))
			{
				cohortKSD <- append(cohortKSD,kstdf$statistic);
				cohortKSP <- append(cohortKSP,kstdf$p.value);
				cohortZKSP <- append(cohortZKSP,kstZdf$p.value);
				cohortZKSD <- append(cohortZKSD,kstZdf$statistic);
			}
			else
			{
				cohortKSD <- append(cohortKSD,NA);
				cohortKSP <- append(cohortKSP,NA);
				cohortZKSP <- append(cohortZKSP,NA);
				cohortZKSD <- append(cohortZKSD,NA);
			}
			
			if (!is.na(kendcor[[1]])) 
			{
				kendall.r <- append(kendall.r,kendcor$estimate);
				kendall.p <- append(kendall.p,kendcor$p.value);
				pearson.r <- append(pearson.r,pearcor$estimate);
				spearman.r <- append(spearman.r,speacor$estimate);
				cStatCorr <- append(cStatCorr,cstat[1]);
			}
			else
			{
				kendall.r <- append(kendall.r,NA);
				kendall.p <- append(kendall.p,NA);
				pearson.r <- append(pearson.r,NA);
				spearman.r <- append(spearman.r,NA);
				cStatCorr <- append(cStatCorr,NA);
			}
			
			
			if (is.null(sizecoef) || is.na(lmodel$coef[sizecoef])) 
			{
				test=NA;
				if (uniType=="Binary")
				{
					IDI <- append(IDI,test);
					NRI <- append(NRI,test);
					zIDI <- append(zIDI,test);
					zNRI <- append(zNRI,test);
					ROCAUC <- append(ROCAUC,test);
					caseN_Z_Low_Tail <- append(caseN_Z_Low_Tail,test);
					caseN_Z_Hi_Tail <- append(caseN_Z_Hi_Tail,test);
					controlN_Z_Low_Tail <- append(controlN_Z_Low_Tail,test);
					controlN_Z_Hi_Tail <- append(controlN_Z_Hi_Tail,test);
				}
				ZGLM <- append(ZGLM,test);
				NeRI <- append(NeRI,test);
				BinRes.p <- append(BinRes.p,test);
				WilcoxRes.p <- append(WilcoxRes.p,test);
				TstudentRes.p <- append(TstudentRes.p,test);
				FRes.p <- append(FRes.p,test);
			}
			else
			{
				if (FullAnalysis)
				{
					if (uniType=="Binary")
					{
						spredict <- predictForFresa(lmodel,data, 'prob');
#						iprob <- improveProb(basepredict,spredict,data[,Outcome]);
						iprob <- .Call("improveProbCpp",basepredict,spredict,data[,Outcome],0);
						IDI <- append(IDI,iprob$idi);
						NRI <- append(NRI,iprob$nri);
						zIDI <- append(zIDI,iprob$z.idi);
						zNRI <- append(zNRI,iprob$z.nri);
						if (length(data[,Outcome])==length(spredict))
						{
							ROCAUC <- append(ROCAUC,pROC::roc( data[,Outcome], spredict,plot=FALSE,auc=TRUE)$auc[1]);
						}
						else 
						{
							ROCAUC <- append(ROCAUC,NA);
						}
						caseN_Z_Low_Tail <- append(caseN_Z_Low_Tail,caseCount1);
						caseN_Z_Hi_Tail <- append(caseN_Z_Hi_Tail,caseCount2);
						controlN_Z_Low_Tail <- append(controlN_Z_Low_Tail,controlCount1);
						controlN_Z_Hi_Tail <- append(controlN_Z_Hi_Tail,controlCount2);
					}

					varResiduals <- residualForFRESA(lmodel,data,Outcome);
#					rprob <- improvedResiduals(baseResiduals,varResiduals);
					rprob <- .Call("improvedResidualsCpp",baseResiduals,varResiduals," ",0);
					NeRI <- append(NeRI,rprob$NeRI);
					BinRes.p <- append(BinRes.p,rprob$p.value);
					WilcoxRes.p <- append(WilcoxRes.p,rprob$WilcoxP.value);
					TstudentRes.p <- append(TstudentRes.p,rprob$tP.value);
					FRes.p <- append(FRes.p,rprob$FP.value);
				}
				else
				{
					if (uniType=="Binary")
					{
						spredict <- predictForFresa(lmodel,data, 'prob');
						iprob <- .Call("improveProbCpp",basepredict,spredict,data[,Outcome],0);
						zIDI <- append(zIDI,iprob$z.idi);
					}
				}
				ZGLM  <- append(ZGLM,abs(modcoef[sizecoef,zcol]));
			}
		}
	}
	
	if (FullAnalysis)
	{
		if (uniType=="Binary")
		{
			orderframe <- data.frame(Name,parent,descrip,cohortMean,cohortStd,cohortKSD,cohortKSP,caseMean,
			caseStd,caseKSD,caseKSP,caseZKSD,caseZKSP,controlMean,controlStd,controlKSD,controlKSP,controlZKSD,
			controlZKSP,t.Rawvalue,t.Zvalue,wilcox.Zvalue,ZGLM,zNRI,zIDI,ROCAUC,cStatCorr,NRI,IDI,NeRI,kendall.r,
			kendall.p,BinRes.p,TstudentRes.p,WilcoxRes.p,FRes.p,caseN_Z_Low_Tail,caseN_Z_Hi_Tail,controlN_Z_Low_Tail,controlN_Z_Hi_Tail);
			switch(rankingTest,
				zIDI=
				{
					orderframe <- with(orderframe,orderframe[order(-zIDI),]);
				},
				zNRI=
				{
					orderframe <- with(orderframe,orderframe[order(-zNRI),]);
				},
				IDI=
				{
					orderframe <- with(orderframe,orderframe[order(-IDI),]);
				},
				NRI=
				{
					orderframe <- with(orderframe,orderframe[order(-NRI),]);
				},
				NeRI=
				{
					orderframe <- with(orderframe,orderframe[order(-NeRI),]);
				},
				Ztest=
				{
					orderframe <- with(orderframe,orderframe[order(-ZGLM),]);
				},
				AUC=
				{
					orderframe <- with(orderframe,orderframe[order(-ROCAUC),]);
				},
				Kendall=
				{
					orderframe <- with(orderframe,orderframe[order(kendall.p),]);
				},
				{
					orderframe <- with(orderframe,orderframe[order(-ROCAUC),]);
				}
			)

		}
		else
		{
			orderframe <- data.frame(Name,parent,descrip,cohortMean,cohortStd,cohortKSD,cohortKSP,cohortZKSD,cohortZKSP,ZGLM,
			NeRI,cStatCorr,spearman.r,pearson.r,kendall.r,kendall.p,BinRes.p,TstudentRes.p,WilcoxRes.p,FRes.p);
			switch(rankingTest,
				NeRI=
				{
					orderframe <- with(orderframe,orderframe[order(-NeRI),]);
				},
				Ztest=
				{
					orderframe <- with(orderframe,orderframe[order(-ZGLM),]);
				},
				CStat=
				{
					orderframe <- with(orderframe,orderframe[order(-cStatCorr),]);
				},
				Kendall=
				{
					orderframe <- with(orderframe,orderframe[order(kendall.p),]);
				},
				{
					orderframe <- with(orderframe,orderframe[order(-cStatCorr),]);
				}
			)
		}

	}
	else
	{
		if (uniType=="Binary")
		{
			orderframe <- data.frame(Name,parent,descrip,ZGLM,zIDI);
			switch(rankingTest,
				zIDI=
				{
					orderframe <- with(orderframe,orderframe[order(-zIDI),]);
				},
				Ztest=
				{
					orderframe <- with(orderframe,orderframe[order(-ZGLM),]);
				}
			)
		}
		else
		{
			orderframe <- data.frame(Name,parent,descrip,ZGLM);
			orderframe <- with(orderframe,orderframe[order(-ZGLM),]);
		}
	}

	row.names(orderframe) <- orderframe$Name;

	# if ((categorizationType == "ERaw")&&((rankingTest=="zIDI")||(rankingTest=="Ztest")))
	# {
		  # isna = is.na(orderframe[,rankingTest])
		  # totb <- sum(orderframe[!isna,rankingTest] > abs(qnorm(0.05/ncol(data))))
		  # varlist <- orderframe[1:totb,1:2]
		  # c.var <- listTopCorrelatedVariables(varlist,data,pvalue=0.05/ncol(data),corthreshold=0.90,method ="spearman")
		  # if (nrow(c.var$correlated.variables)>0)
		  # {
			# addednames <- as.data.frame(c.var$correlated.variables)
			# rownames(orderframe) <- orderframe$Name
			# lastorder <- orderframe[rownames(addednames),]
			# lastorder$Name <- addednames[,1]
			# lastorder$parent <- rownames(addednames)
			# rownames(lastorder) <- paste("D_",rownames(addednames),sep = "")
			# colnames(lastorder) <- colnames(orderframe)
			# orderframe <- rbind(orderframe,lastorder)
			# switch(rankingTest,
				   # zIDI=
					# {
					  # orderframe <- with(orderframe,orderframe[order(-zIDI),]);
					# },
					# Ztest=
					# {
					  # orderframe <- with(orderframe,orderframe[order(-ZGLM),]);
					# }
				# )
		  # }
	# }
	
	
	result <- list(orderframe=orderframe,
	variableList=variableList,
	formula=formula,
	Outcome=Outcome,
	data=data,
	categorizationType=categorizationType,
	type=type,
	rankingTest=rankingTest,
	cateGroups=cateGroups,
	raw.dataFrame=raw.dataFrame,
	description=description,
	uniType=uniType)
	

	return (result);
	
}
