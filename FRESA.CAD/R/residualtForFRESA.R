residualForFRESA <-
function (object,testData,Outcome,eta=0.05) 
{
#Set eta to zero to get prediction residuals for cox models. eta to 1 get the Martingale residuals
	classlen=length(class(object))
	
	cobj <- substr(class(object)[classlen], 1, 2);
	if ((cobj!="co") && !is.null(object$family[1])) cobj = "glm";
	switch(cobj,
		co =
		{
			if (classlen==1)
			{
				lpp <- 1.0/(1.0+exp(-predict(object,testData,'lp',na.action=na.omit)));
				s <- is.na(lpp);
				if (any(s)) 
				{
					lpp[s] = 1.0e10; # set to a large residual
				}
				lppres <- lpp - testData[,Outcome];
				matingale <- testData[,Outcome]-predict(object,testData,'expected',na.action=na.omit);
				s <- is.na(matingale);
				if (any(s))
				{
					matingale[s] = 1.0e10; # set to a large residual
				}
				out <- (1-eta)*lppres - eta*matingale;			
			}
			else
			{
				out <- (1 - 2*testData[,Outcome]);
			}
		},
		lm =
		{
			out <- predictForFresa(object,testData,'linear') - testData[,Outcome];
		},
		fi =
		{
			if (object$type == "LM")
			{
				out <- predictForFresa(object,testData,'linear') - testData[,Outcome];
			}
			else
			{
				out <- predictForFresa(object,testData,'prob') - testData[,Outcome];
			}
		},
		{
#			if (object$family[1] == "binomial")
#			{
#				out <- predictForFresa(object,testData,'prob') - testData[,Outcome];
#			}
#			else
#			{
#				out <- predictForFresa(object,testData,'linear') - testData[,Outcome];
#			}
#		out <- predict(object,testData,'response',na.action=na.omit) - testData[,Outcome];
		out <- 0 - testData[,Outcome];
		}
	)	
	s <- is.na(out);
	if (any(s)) 
	{
		cat("Warning NA predictFor NeRIs \n");
		out[s] <- 1.0e10; # Set a large residual
	}
    return (out)
}
