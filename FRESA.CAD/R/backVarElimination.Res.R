backVarElimination_Res <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),setIntersect=1,adjsize=1) 
{

	back.var.NeRISelection <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),setIntersect=1) 
	{
		type <- match.arg(type);
	  
		varsList <- as.list(attr(terms(object),"variables"))
		termList <- as.list(attr(terms(object),"term.labels"))
		
		modsize <- length(termList);
		
		removeID = 0;

		outCome = paste(varsList[2]," ~ ");
		if (setIntersect==0) 
		{
			outCome = paste(outCome," 0  ");
		}
		else
		{
			outCome = paste(outCome," 1  ");
		}
		startlist = 3 ;
		frm1 = outCome;
		cpv = pvalue; 
		if (length(varsList)>=startlist)
		{
			for ( i in startlist:length(varsList))
			{
				frm1 <- paste(frm1,paste(" + ",varsList[i]));
			}
		}
#		cat ("Len: ",length(termList)," : ",frm1,"\n")
		ftmp <- formula(frm1);
		bckform <- frm1;
		startSearch = startlist + startOffset;
		if (length(termList)>1)
		{
			FullModel <- modelFitting(ftmp,data,type)
			FullResiduals <- residualForFRESA(FullModel,data,Outcome);

			if (length(varsList)>startSearch)
			{
				for ( i in startSearch:length(varsList))
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
					redModel <- modelFitting(ftmp,data,type,TRUE)
					if (inherits(redModel, "try-error"))
					{
						redModel <- FullModel
					}
					

					redResiduals <- residualForFRESA(redModel,data,Outcome);
					iprob <- improvedResiduals(redResiduals,FullResiduals,testType);
					if (iprob$p.value>cpv)
					{
						cpv = iprob$p.value;
						removeID = i;
					}
				}
			}
			if ((length(varsList) == startSearch) && (removeID == startSearch)) 
			{
				removeID = -1;
			}

			if (removeID>0)
			{
				frm1 = outCome;
				for ( i in startlist:length(varsList))
				{
					if (i != removeID)
					{
						frm1 = paste(frm1,paste(" + ",varsList[i]));
					}
				}
				ftmp <- formula(frm1);
				bckform <- frm1;
#				cat("formula :",frm1,"\n")
				FullModel <- modelFitting(ftmp,data,type)
			}
		}
		else
		{
			FullModel <- object;
		}
#		cat("Reduced: \n")
#		print(summary(FullModel));
		result <- list(Model=FullModel,Removed=removeID,backfrm=bckform);

		return (result)
	}

	bkobj <- NULL;
	if (adjsize>1)
	{
		bkobj <- backVarElimination_Res(object,pvalue,Outcome,data,startOffset,type,testType,setIntersect,adjsize=1); # remove features that do not improve residuals
		object <- bkobj$back.model;
		adjsize = floor(adjsize);
		adjsize <- min(adjsize,ncol(data)-1);
#		cat("Adjusted Size:",adjsize,"\n");
	}

	changes=1;
	loops=0;
    model <- object;
	beforeFSCmodel <- object;
	changes2<-0
	while ((changes>0) && (loops<100))
	{
		p.elimin <- pvalue;
		if (adjsize>0)
		{		
			modsize <- length(as.list(attr(terms(model),"term.labels")));	
			if (modsize<1) modsize=1;
			qvalue <- 4*pvalue;
			if (qvalue < 0.1) qvalue=0.1 # lests keep a the minimum q-value to 0.1
			p.elimin <- min(pvalue,modsize*qvalue/adjsize) # BH alpha  the elimination p-value
		}

		bk <- back.var.NeRISelection(model,p.elimin,Outcome=Outcome,data=data,startOffset,type,testType,setIntersect);
#		cat("Used p :",p.elimin,"Formula <- ", bk$backfrm,"\n");
		changes = as.integer(bk$Removed);
		if (changes>0)
		{
		  loops = loops + 1;
			changes2<- as.character(as.list(attr(terms(model),"variables")))[which(!(as.character(as.list(attr(terms(model),"variables")))%in%as.character(as.list(attr(terms(bk$Model),"variables")))))]
			model = bk$Model;
			if (length(changes2)>1)
			{
				changes2<-changes2[2]
			}
				if (adjsize>1)
				{
					weight <- 1.0;
					if ((length(beforeFSCmodel$coefficients)>0)&&(length(model$coefficients)>0))
					{
						for (i in 1:length(beforeFSCmodel$coefficients))
						{
							notadded = TRUE;
							for (j in 1:length(model$coefficients))
							{
								if (names(beforeFSCmodel$coefficients)[i] == names(model$coefficients)[j])
								{
									beforeFSCmodel$coefficients[i] <- (weight*beforeFSCmodel$coefficients[i] + (1-weight)*model$coefficients[j]); # it will average the two
									notadded=FALSE;
								}
							}
							if (notadded)
							{
								beforeFSCmodel$coefficients[i] <- weight*beforeFSCmodel$coefficients[i]; # it will average with zero
							}
						}
					}
				}
		}
		if (changes < 0)
		{
			changes2<- changes
		}
		model = bk$Model;
	}
#	cat("Reduced Model:",bk$backfrm,"\n")
#	print(summary(model));

	modelReclas <- getVar.Res(model,data=data,Outcome=Outcome,type=type);
	
	result <- list(back.model= model,
	loops=loops,
	reclas.info=modelReclas,
	back.formula=formula(bk$backfrm),
	bootCV=NULL,
	lastRemoved=changes2,
	number.of.independent=adjsize,
	beforeFSC.model=beforeFSCmodel,
	beforeFSC.formula=formula(bkobj$backfrm));
	return (result);
}