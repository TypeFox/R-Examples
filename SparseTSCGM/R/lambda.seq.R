
lambda.seq <- function(SS, SA,nlambda)
   {
	  if(is.null(nlambda)) nlambda=10
	  else if(!is.null(nlambda)) nlambda=nlambda
    lambda.min.ratio=0.1
	  d= dim(SS)[2]
		lambda.max1 = max(max(SS-diag(d)),-min(SS-diag(d)))
		lambda.min1 = lambda.min.ratio*lambda.max1
		lambda1 = exp(seq(log(lambda.max1), log(lambda.min1), length = nlambda))
		 lambda.min.ratio2=0.15
		lambda.max2 = max(max(SA),-min(SA))
		lambda.min2 = lambda.min.ratio2*lambda.max2
		lambda2 = exp(seq(log(lambda.max2), log(lambda.min2), length = nlambda))
		return(list(lambda1=lambda1, lambda2=lambda2))
	}

