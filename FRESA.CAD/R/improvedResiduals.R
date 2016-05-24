improvedResiduals <-
function (oldResiduals,newResiduals,testType=c("Binomial","Wilcox","tStudent","Ftest"))
{
	testType <- match.arg(testType)

# initialize the variables
	pwil <- 1;
	pbin <- 1;
	ptstu <- 1;
	f.test <- 1;
	p1=0;
	p2=0;
	improved = 0;
	
	
	# compute the difference in residuals 
	size <- length(oldResiduals);
	size1 <- size - 1;

	oldres <- abs(oldResiduals);
	newres <- abs(newResiduals);

	delta <- newres - oldres;

	# Count the number of samples that improved the residuals 
	reduction <- sum(delta < 0);

	# Count the number of samples that worsen the residuals 
	increase  <- sum(delta > 0);

	improved <- (reduction-increase)/size; 												# the net improvement in residuals
	p1 <- reduction/size;																#proportion of subjects with improved residuals
	p2 <- increase/size;																#proportion of subjects with worst residuals
	rss1 <- sum(oldResiduals^2)
	rss2 <- sum(newResiduals^2)/size1;
	f.test <- 1-pf(rss1/rss2-size1,1,size1); 											# Just compare an improvement of the residual Variance
	tdat <- try(t.test(oldres, newres, paired = TRUE,alternative = "greater"))
	if (!inherits(tdat, "try-error"))
	{
		ptstu <- tdat$p.value  	# let compute that the probability that the new residuals are better than the old residuals via t-test
	}
	else 
	{
		ptstu <- 1.0;
	}
	if (improved>0)
	{
		pbin <- binom.test(reduction,size,alternative = "greater")$p.value 					# Lets do a sign test to test a significant improvement in residual variance
		pwil <- wilcox.test(oldres, newres, paired = TRUE,alternative = "greater")$p.value  # let compute that the probability that the new residuals are better than the old residuals via wilcoxon
	}
	switch(testType, 
	    Binomial =
        {
			pvalue = pbin;
		},
		Wilcox =
		{
			pvalue = pwil;
		},
		tStudent =
		{
			pvalue = ptstu;
		},
		Ftest =
		{
			pvalue = f.test;
		}
	)
		

	result <- list(p1=p1,
	p2=p2,
	NeRI=improved,
	p.value = pvalue,
	BinP.value = pbin,
	WilcoxP.value = pwil,
	tP.value = ptstu,
	FP.value = f.test);
	return (result)
}