CVonePair <-function(X,y,nFolds,foldId,hyperpara=c(1,0.1),Epis=FALSE, 
					prior=c("lassoNEG","lasso","elastic net"), family = c("gaussian","binomial"), verbose = 0, group = FALSE)
{	
	if(prior=="lassoNEG")
	{
		pr =1;
	}else if(prior=="lasso")
	{
		pr = 2;
	}else
	{
		pr =3;
	}
	if(family =="gaussian")
	{
		model =0;
	}else
	{
		model = 1;
	}
	N = nrow(X);
	K = ncol(X);
	nLogL = rep(0,4);
	output<-.C("cvOnePara",
			BASIS 		= as.double(X),
			y 		= as.double(y),
			foldId  	= as.integer(foldId),
			nfolds 		= as.integer(nFolds),
			n  	= as.integer(N),
			k 		= as.integer(K),
			verbose =as.integer(verbose),
			hyperpara 		= as.double(hyperpara),
			nLogL  	= as.double(nLogL),
			epistasis 		= as.integer(Epis),
			pr  	= as.integer(pr),			
			glm	= as.integer(model),
			group = as.integer(group),
			PACKAGE 	="EBglmnet");
output$nLogL #negative log likelihood
}

