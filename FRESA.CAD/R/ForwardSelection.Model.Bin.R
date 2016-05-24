ForwardSelection.Model.Bin <-
function(size=100,fraction=1.0,pvalue=0.05,loops=100,covariates="1",Outcome,variableList,data,maxTrainModelSize=10,type=c("LM","LOGIT","COX"),timeOutcome="Time",selectionType=c("zIDI","zNRI","Both"),loop.threshold=20,interaction=1,cores=4)
{
#	    R_CStackLimit = -1;
	    type <- match.arg(type)
	    seltype <- match.arg(selectionType)


	  	Outcome<-as.character(Outcome);

	  	if (type=="COX")
	  		timeOutcome<-as.character(timeOutcome)
	  	else
	  		timeOutcome="";
	    if (is.na(size))
	    {
	      stop("Number of variables to be used is not defined\n")
	    }

		acovariates <- covariates[1];
		if (length(covariates)>1)
		{
			for (i in 2:length(covariates))
			{	
				acovariates <- paste(acovariates," + ",covariates[i])
			}
		}

		
		vnames <- as.vector(variableList[,1]);
		
		if (nrow(variableList)>1)
		{
			if (nrow(variableList)<size) size = nrow(variableList)
			frm <- paste(Outcome,"~",acovariates," + ",timeOutcome);
			for (i in 1:size)
			{
				frm <- paste(frm," + ",vnames[i])
			}
#			cat(frm,"\n")
			modelFrame <- model.frame(formula(frm),data);
		}
		else
		{
			modelFrame <- data;	
		}

#		cat(frm,"\n")
		
	    colNames=colnames(modelFrame);
	    
	    output<-.Call("ReclassificationFRESAModelCpp",size, fraction, pvalue, loops, covariates, Outcome,as.vector(variableList[,1]), maxTrainModelSize, type, timeOutcome, seltype,loop.threshold, interaction,data.matrix(modelFrame),colNames,cores);

		rloops = min(loops,100);
		rpvalue = min(pvalue,0.01);
	    if (loops>1) 
		{
			randoutput <-.Call("ReclassificationFRESAModelCpp",size, fraction, rpvalue, rloops, covariates, "RANDOM" ,as.vector(variableList[,1]), maxTrainModelSize, type, timeOutcome, seltype,loop.threshold, interaction,data.matrix(modelFrame),colNames,cores);
		}
		else 
		{
			randoutput <-output;
			rpvalue=pvalue;
		}
		
		zthr = abs(qnorm(pvalue)); 
		zthrO = abs(qnorm(pvalue*pvalue));
		zthr2 = zthr;
		if (fraction<1) 
		{
			zthr2 = zthr*sqrt(fraction);
		}
		if (zthr2<abs(qnorm(0.1))) zthr2 = abs(qnorm(0.1));
	


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
#			cat(output$formula.list[i],"\n")
			avgsize = avgsize+ str_count(randoutput$formula.list[i],"\\+") - 1;
		}
		avgsize = (pvalue/rpvalue)*(nrow(variableList)/size)*(avgsize/rloops);
		cat ("To Test Variables:",nrow(variableList),"# Variables:",size," Average size =",avgsize,"\n");
		if (loops==1) avgsize=0;

			mynames <- output$mynames + 1 
			topvar <- table(mynames);
	#		print(mynames)
			if (length(topvar)>1)
			{
				topvar <- topvar[order(-topvar)];
			}
			topvarID <- as.numeric(rownames(topvar));
	#		print(topvar)


			
			frm1 <- baseForm;
			frm1 <- paste(frm1," + ");
			frm1 <- paste(frm1,vnames[topvarID[1]]);
			ftmp <- formula(frm1);
	#		cat(frm1," <- Start Formula \n")
			bestmodel <- modelFitting(ftmp,data,type)
			
			if ( !inherits(bestmodel, "try-error"))
			{
				bestpredict <- predictForFresa(bestmodel,data,'prob');

				vnames_model <- vector();
				model_zmin <- vector();
				varlist <- vector();
				inserted = 1;
				
				vnames_model <- append(vnames_model,vnames[topvarID[1]]);
				varlist <- append(varlist,topvarID[1]);
				model_zmin <- append(model_zmin,NA);
				
				if (length(topvar)>1)
				{
					for ( i in 2:length(topvar))
					{

						if (loops > loop.threshold) 
						{
							frec <- topvar[i];
							if (!is.na(frec))
							{
								if ((frec/loops) < 1.0/(2*loop.threshold+1.0)) 
								{ 
									topvar[i] = 0;
								}
							}
						}

						if ((topvar[i] > 0) && (inserted < maxTrainModelSize))
						{
							kinserted = 0
							kins = 0 
							frma <- paste(frm1," + ");
							frma <- paste(frma,vnames[topvarID[i]]);
	#						cat(frma,"\n");
							ftmp <- formula(frma);
							newmodel <- modelFitting(ftmp,data,type,TRUE)
							if ( !inherits(newmodel, "try-error"))
							{
	#							iprob <- improveProb(bestpredict,predictForFresa(newmodel,data,'prob'),data[,Outcome]);
								iprob <- .Call("improveProbCpp",bestpredict,predictForFresa(newmodel,data,'prob'),data[,Outcome],0);
								if (seltype=="zIDI") 
								{
									zmin = iprob$z.idi;
								}
								else
								{
									zmin = iprob$z.nri;
								}
								if (is.numeric(zmin) && !is.na(zmin) && (zmin>zthr))
								{
									bestpredict <- predictForFresa(newmodel,data,'prob');
									frm1 <- frma;
									vnames_model <- append(vnames_model,vnames[topvarID[i]]);
									model_zmin <- append(model_zmin,zmin);
									varlist <- append(varlist,topvarID[i]);
									inserted = inserted + 1;
									kins = 1;
								}	
								if (interaction == 2)
								{
									for (nlist in 1:inserted)
									{
										if (kins==1)
										{
											zthrOl=zthr;
											frma <- paste(frm1," + I(",vnames[varlist[nlist]],"*",vnames[topvarID[i]],")")
										}
										else
										{
											frma <- paste(frm1," + ",vnames[topvarID[i]]," + I(",vnames[varlist[nlist]],"*",vnames[topvarID[i]],")")
											zthrOl=zthrO;
										}
										ftmp <- formula(frma);
										newmodel <- modelFitting(ftmp,data,type,TRUE)
										if ( !inherits(newmodel, "try-error"))
										{
	#										iprob <- improveProb(bestpredict,predictForFresa(newmodel,data,'prob'),data[,Outcome]);
											iprob <- .Call("improveProbCpp",bestpredict,predictForFresa(newmodel,data,'prob'),data[,Outcome],0);
											if (seltype=="zIDI") 
											{
												zmin = iprob$z.idi;
											}
											else
											{
												zmin = iprob$z.nri;
											}
											if (is.numeric(zmin) && !is.na(zmin) && (zmin>zthrOl))
											{
												bestpredict <- predictForFresa(newmodel,data,'prob');
												frm1 <- frma;
												vnames_model <- append(vnames_model,vnames[topvarID[i]]);
												model_zmin <- append(model_zmin,zmin);
												kinserted = kinserted + 1;
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
				}
			}
	#		barplot(topvar);
	#		titname <- paste ( "Var Frequency Completed");
	#		title(main=titname);
			ftmp <- formula(frm1);
			bestmodel <- modelFitting(ftmp,data,type)

		#	print(summary(bestmodel));
		

		result <- list(final.model=bestmodel,
		var.names=vnames_model,
		formula=ftmp,
		ranked.var=topvar,
		z.selectionType=model_zmin,
		formula.list=output$formula.list,
		average.formula.size=avgsize,
		variableList=variableList
		);
#		cat ("Final :",frm1,"\n")
	    return (result);
}
	