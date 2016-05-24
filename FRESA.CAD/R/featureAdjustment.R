featureAdjustment <-
function(variableList,baseModel,strata=NA,data,referenceframe,type=c("LM","GLS"),pvalue=0.05,correlationGroup = "ID") 
{

if (!requireNamespace("nlme", quietly = TRUE)) {
   install.packages("nlme", dependencies = TRUE)
} 


	type <- match.arg(type);
	## the reference frame will be used to predict a variable from the basemodel. At output the residuals are returned.
	## strata is a numeric column varname in the data frame from 0 to S, where S is the maximum number of strata
	colnamesList <- as.vector(variableList[,1]);
	size = length(colnamesList);
	if (!is.na(strata)) 
	{
		maxStrata = max(referenceframe[,strata],na.rm = TRUE);
		minStrata = min(referenceframe[,strata],na.rm = TRUE);
	}
	else 
	{
		maxStrata=1;
		minStrata=1;
	}
	cat ("Min Strata:",minStrata,"Max Strata:",maxStrata,"\n");

#	SortedCtr = referenceframe;
#	nrowsCtr =  nrow(referenceframe);
	created = 0;
	for (sta in minStrata:maxStrata)
	{
		if (!is.na(strata))
		{
			stracondition = paste (strata,paste('==',sta));
			strastatement = paste ("subset(referenceframe,",paste(stracondition,")"));
			cat ("Strata:",stracondition,"\n");
			cstrataref <- eval(parse(text=strastatement));
			strastatement = paste ("subset(data,",paste(stracondition,")"));
			cstrata <- eval(parse(text=strastatement));
			cat ("Rows:",nrow(cstrataref),"Rows 2",nrow(cstrata)," \n");
		}
		else
		{
			cstrataref = referenceframe;
			cstrata = data;
		}
		if ((nrow(cstrata)>1) && ( nrow(cstrataref)>1))
		{
			for (i in 1:size)
			{ 
				ftm1 <- paste(colnamesList[i],paste(" ~ ",baseModel));
				ftmp <- formula(ftm1);

				switch(type, 
					LM = 
					{ 
						model <- lm(ftmp,data=cstrataref,na.action=na.exclude)
						f <- summary(model)$fstatistic
						f1 = f[1];
						p <- pf(f[1],f[2],f[3],lower.tail=FALSE)			
					},
					GLS =
					{
						model <- eval(parse(text=paste("try(nlme::gls(formula(",ftm1,"),cstrataref,na.action=na.exclude,correlation = nlme::corAR1(form = ~ 1 | ",correlationGroup,")))")))
						dgf = nrow(cstrataref)-1;
						rss1 <- var(cstrataref[,colnamesList[i]],na.rm = TRUE)
						rss2 <- var(model$residuals,na.rm = TRUE);
						f1 = rss1/rss2;
						p <- 1-pf(dgf*rss1/rss2-dgf,1,dgf); 
					},
					{ 
						model <- lm(ftmp,data=cstrataref,na.action=na.exclude)
						f <- summary(model)$fstatistic
						f1 = f[1];
						p <- pf(f[1],f[2],f[3],lower.tail=FALSE)			
					}
				)
		
				cat(" Variable:\t ",colnamesList[i],"\t F Stats:\t ",f1,"\t P-value:\t",p,"\n");
				if (!is.na(p))
				{
					switch(type, 
						LM = 
						{ 
							if (p<pvalue)
							{
								cstrata[,colnamesList[i]] <- cstrata[,colnamesList[i]]-predict(model,cstrata);
							}
							else
							{
								avg <- mean(cstrataref[,colnamesList[i]],na.rm = TRUE);
								cstrata[,colnamesList[i]] <- cstrata[,colnamesList[i]]-avg;
							}
						},
						GLS =
						{
							if (p<pvalue)
							{
								avg <-model$coef[1];
								cstrata[,colnamesList[i]] <- cstrata[,colnamesList[i]]-predict(model,cstrata)+avg;
							}
						}
					)
				}
			}
			if (created == 1) 
			{
				AdjustedFrame <- rbind(AdjustedFrame,cstrata);
			}
			else
			{
				created = 1;
				AdjustedFrame = cstrata;
			}
		}
	}
	for (i in 1:size)		
	{ 
		var1 <- var(data[,colnamesList[i]],na.rm = TRUE);
		var2 <- var(AdjustedFrame[,colnamesList[i]],na.rm = TRUE);
		cat(" Variable: \t",colnamesList[i],"\t Var Ini: \t",var1,"\t Var End:\t",var2,"\t F:\t",var1/var2,"\n");
	}
	if (created ==0 ) 
	{
		AdjustedFrame=NULL;
	}

	return (AdjustedFrame);
}
