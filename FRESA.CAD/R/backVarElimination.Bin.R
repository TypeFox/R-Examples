backVarElimination_Bin <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),selectionType=c("zIDI","zNRI"),adjsize=1) 
{
  	seltype <- match.arg(selectionType)

	back.var.IDISelection <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),selectionType=c("zIDI","zNRI")) 
	{
		seltype <- match.arg(selectionType)
		type <- match.arg(type);
	  
		varsList <- as.list(attr(terms(object),"variables"))
		
		
		cthr = abs(qnorm(pvalue));
		removeID = 0;

		outCome = paste(varsList[2]," ~ ");
		startlist = 3 ;
		frm1 = outCome;
		if (length(varsList)>=startlist)
		{
			for ( i in startlist:length(varsList))
			{
				frm1 <- paste(frm1,paste(" + ",varsList[i]));
			}
		}
		else
		{
			frm1 <- paste(frm1," 1 ")
		}
		
		ftmp <- formula(frm1);
		bckform <- frm1;
		FullModel <- modelFitting(ftmp,data,type)
		startSearch = startlist + startOffset;
		if ( !inherits(FullModel, "try-error"))
		{
			FullPredict <- predictForFresa(FullModel,data,'prob');
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
					if ( !inherits(redModel, "try-error"))
					{
						redPredict <- predictForFresa(redModel,data,'prob');
#						iprob <- improveProb(redPredict,FullPredict,data[,Outcome]);
						iprob <- .Call("improveProbCpp",redPredict,FullPredict,data[,Outcome],0);
						if (seltype=="zIDI") 
						{
							ztst = iprob$z.idi;
						}
						else
						{
							ztst = iprob$z.nri;
						}
						if (is.na(ztst)) ztst=0;
						if (ztst<cthr)
						{
							cthr = ztst;
							removeID = i;
						}
					}
				}
			}
		}

		if ((length(varsList) == startSearch) && (removeID == startSearch)) 
		{
			removeID = -1;
		}
		
		if (removeID > 0)
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
			FullModel <- modelFitting(ftmp,data,type)
		}
#		cat("removed : ",removeID,"Final Model: \n")
#		print(summary(FullModel));
		result <- list(Model=FullModel,Removed=removeID,backfrm=bckform);

		return (result)
	}

	bkobj <- NULL;
	if (adjsize>1)
	{
		bkobj <- backVarElimination_Bin(object,pvalue,Outcome,data,startOffset,type,selectionType,adjsize=1); 
		object <- bkobj$back.model;
		adjsize = floor(adjsize);
		adjsize <- min(adjsize,ncol(data)-1);
#		cat("Adjusted Size:",adjsize,"\n");
	}

	changes=1;
	loops=0;
    model <- object;
	beforeFSCmodel <- object;
	mydataFrame <- data;
	myOutcome <- Outcome;
	changes2 <- 0
	while ((changes>0) && (loops<100))
	{
		p.elimin <- pvalue;
		if (adjsize>1)
		{
			modsize <- length(as.list(attr(terms(model),"term.labels")));	
			if (modsize<1) modsize=1;
			qvalue <- 4*pvalue;
			if (qvalue < 0.1) qvalue=0.1 # lests keep a the minimum q-value to 0.1
			p.elimin <- min(pvalue,modsize*qvalue/adjsize) # BH alpha the elimination p-value
		}

		bk <- back.var.IDISelection(model,p.elimin,Outcome=myOutcome,data=mydataFrame,startOffset,type,seltype);

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
#	print(summary(model));
#	cat("Reduced Model:",bk$backfrm,"\n")
	modelReclas <- getVar.Bin(model,data=mydataFrame,Outcome=myOutcome,type);
	result <- list(back.model=model,
	loops=loops,
	reclas.info=modelReclas,
	back.formula=formula(bk$backfrm),
	lastRemoved=changes2,
	beforeFSC.model=beforeFSCmodel,
	beforeFSC.formula=formula(bkobj$backfrm));
	return (result);
}