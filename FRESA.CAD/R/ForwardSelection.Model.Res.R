ForwardSelection.Model.Res <-
function(size=100,fraction=1.0,pvalue=0.05,loops=100,covariates="1",Outcome,variableList,data,maxTrainModelSize=10,type=c("LM","LOGIT","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),timeOutcome="Time",loop.threshold=20,interaction = 1,cores = 4)
{
#	R_CStackLimit = -1;

	if (is.na(size))
	{
		stop("Size: Number of variables to be explored is not defined\n")
	}


	type <- match.arg(type)
	testType <- match.arg(testType)
	Outcome<-as.character(Outcome);

	if (type=="COX")
		timeOutcome<-as.character(timeOutcome)
	else
		timeOutcome="";

	vnames <- as.vector(variableList[,1]);
	acovariates <- covariates[1];
	if (length(covariates)>1)
	{
		for (i in 2:length(covariates))
		{	
			acovariates <- paste(acovariates," + ",covariates[i])
		}
	}

	if (nrow(variableList)>1)
	{
		if (nrow(variableList)<size) size = nrow(variableList)
		frm <- paste(Outcome,"~",acovariates," + ",timeOutcome);
		for (i in 1:size)
		{
			frm <- paste(frm," + ",vnames[i])
		}
		modelFrame <- model.frame(formula(frm),data);
	}
	else
	{
		modelFrame <- data;
	}

	colNames=colnames(modelFrame);
	
	output<-.Call("ForwardResidualModelCpp",size, fraction, pvalue, loops, covariates, Outcome,as.vector(variableList[,1]), maxTrainModelSize, type, timeOutcome, testType,loop.threshold, interaction,data.matrix(modelFrame),colNames,cores);
	randoutput <- output;
	rloops = min(loops,100);
	rpvalue = min(pvalue,0.01);
	if (loops>1) 
	{
		randoutput<-.Call("ForwardResidualModelCpp",size, fraction, rpvalue, rloops, covariates, "RANDOM",as.vector(variableList[,1]), maxTrainModelSize, type, timeOutcome, testType,loop.threshold, interaction,data.matrix(modelFrame),colNames,cores);
	}
	else
	{
		randoutput <- output;
		rpvalue=pvalue;
	}

	mynames <- output$mynames + 1;
	formula.list <- output$formula.list

	pthr = pvalue;
	pthrO = pvalue*pvalue;
	baseForm = Outcome;
#For Cox  models 
	if (type == "COX")
	{
	  baseForm = paste("Surv(",timeOutcome);
	  baseForm = paste(baseForm,paste(",",paste(Outcome,")")));
	}

	baseForm = paste(baseForm,paste(" ~ ",acovariates));

	
	avgsize = 0;
	for (i in 1:rloops)
	{
		avgsize = avgsize+ str_count(randoutput$formula.list[i],"\\+") - 1;
	}
	avgsize = (pvalue/rpvalue)*(nrow(variableList)/size)*(avgsize/rloops);
	cat ("Average size =",avgsize,"\n");

		pthr2 = 1-pnorm(sqrt(fraction)*abs(qnorm(pthr)));
		if (pthr2>0.1) pthr2 = 0.1;


		topvar <- table(mynames);
		
		frm1 <- baseForm;
		vnames_model <- vector();
		model_ziri <- vector();
		if (length(topvar)>1)
		{
			topvar <- topvar[order(-topvar)];
			topvarID <- as.numeric(rownames(topvar));

			frm1 <- paste(frm1," + ");
			frm1 <- paste(frm1,vnames[topvarID[1]]);
			
			ftmp <- formula(frm1);
			bestmodel <- modelFitting(ftmp,data,type)
	#		cat(frm1,"b \n")

			bestResiduals <- residualForFRESA(bestmodel,data,Outcome);

			vnames_model <- append(vnames_model,vnames[topvarID[1]]);
			model_ziri <- append(model_ziri,1);
			varlist <- vector();
			varlist <- append(varlist,topvarID[1]);
			inserted = 1
			for ( i in 2:length(topvar))
			{
				if (loops > loop.threshold) 
				{
					frec <- topvar[i];
					if (!is.na(frec))
					{
						if ((frec/loops) < 1.0/(2.0*loop.threshold+1.0)) 
						{ 
							topvar[i] <- 0;
						}
					}
				}
				if(topvar[i] > 0)
				{
					frma <- paste(frm1," + ");
					frma <- paste(frma,vnames[topvarID[i]]);
	#				cat(frma," b \n");

					
					ftmp <- formula(frma);
					newmodel <- modelFitting(ftmp,data,type,TRUE)
					kins = 0
					if ( !inherits(newmodel, "try-error"))
					{
	#					iprob <- improvedResiduals(bestResiduals,residualForFRESA(newmodel,data,Outcome),testType);
						iprob <- .Call("improvedResidualsCpp",bestResiduals,residualForFRESA(newmodel,data,Outcome),testType,0);
	#					cat(frma," c \n");
						piri <- iprob$p.value;
						if (piri<pthr)
						{
							bestResiduals <- residualForFRESA(newmodel,data,Outcome);
	#						cat(frma," d \n");
							frm1 <- paste(frm1," + ");
							frm1 <- paste(frm1,vnames[topvarID[i]]);
	#						cat(frm1," 1a \n");
							varlist <- append(varlist,topvarID[i]);
							vnames_model <- append(vnames_model,vnames[topvarID[i]]);
							model_ziri <- append(model_ziri,abs(qnorm(piri)));
		#					print(summary(newmodel));
							inserted = inserted + 1;
							kins = 1
						}	
						if (interaction == 2)
						{
							for (nlist in 1:inserted)
							{
								if (kins==1)
								{
									pthrOl=pthr;
									frma <- paste(frm1," + I(",vnames[varlist[nlist]],"*",vnames[topvarID[i]],")")
								}
								else
								{
									frma <- paste(frm1," + ",vnames[topvarID[i]]," + I(",vnames[varlist[nlist]],"*",vnames[topvarID[i]],")")
									pthrOl=pthrO;
								}
								ftmp <- formula(frma);
								newmodel <- modelFitting(ftmp,data,type,TRUE)
								if ( !inherits(newmodel, "try-error"))
								{
	#								print(summary(newmodel));
	#								iprob <- improvedResiduals(bestResiduals,residualForFRESA(newmodel,data,Outcome),testType);
									iprob <- .Call("improvedResidualsCpp",bestResiduals,residualForFRESA(newmodel,data,Outcome),testType,0);
									
									piri <- iprob$p.value;
									if (is.numeric(piri) && !is.na(piri) && (piri<pthrOl))
									{
										bestResiduals <- residualForFRESA(newmodel,data,Outcome);
										frm1 <- frma;
										vnames_model <- append(vnames_model,vnames[topvarID[i]]);
										model_ziri <- append(model_ziri,abs(qnorm(piri)));
										if (kins == 0)
										{
											varlist <- append(varlist,topvarID[i]);
											inserted = inserted + 1;
										}
										kins =1
									}
								}							
							}
						}
					}
				}
			}
	#		barplot(topvar);
	#		titname <- paste ( "Var Frequency Completed");
	#		title(main=titname);
	#		print(topvar)
		}


		ftmp <- formula(frm1);
		bestmodel <- modelFitting(ftmp,data,type)

	#	cat(frm1," Final \n");
	
	result <- list(final.model=bestmodel,
	var.names=vnames_model,
	formula=ftmp,
	ranked.var=topvar,
	z.NeRIs=model_ziri,
	formula.list=formula.list,
	average.formula.size=avgsize,
	variableList=variableList);
	
	return (result);
}
