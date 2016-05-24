modelfit <- function(x, y, nfolds, 
	penalty = c("LASSO", "SCAD", "MCP"),
	family=c("gaussian","binomial")) {
  penalty <- match.arg(penalty)
  family <- match.arg(family)
  y <- drop(y)
  y <- as.numeric(y)
  x <- as.matrix(x)
  if(penalty == "LASSO"){					
    cvfit <- cv.glmnet(x=x,y=y,nfolds=nfolds,alpha=1,
	    		type.measure="deviance",maxit=1e6,family=family)
    coefit<-coef(cvfit,s="lambda.min")
	}
  if(penalty == "MCP" || penalty == "SCAD"){
    cvfit<-cv.ncvreg(X=x,y=y,nfolds=nfolds,penalty=penalty,
	    		family=family,max.iter=1e4)
    coefit<-cvfit$fit$beta[,cvfit$min]
	}
  modelfit <- 1-(coefit[-1]==0)
  list(coefit = coefit, modelfit = modelfit)
}