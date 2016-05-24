modelFitting <-
function (model.formula,data,type = c("LOGIT", "LM","COX"),fast=FALSE,...) 
{

	type <- match.arg(type);
	
#if (!requireNamespace("speedglm", quietly = TRUE)) {
#   install.packages("speedglm", dependencies = TRUE)
#} 

	if (!fast)
	{

		switch(type, 
			LM = 
			{ 
			  fittedModel <- try(lm(model.formula,data=data,na.action=na.exclude));
			 },
			LOGIT =
			{
	          fittedModel <- try(glm(model.formula,data=data,na.action=na.exclude,family=binomial(link=logit),...));
#			  fittedModel <- try(speedglm::speedglm(model.formula,data=data,na.action=na.exclude,family=binomial(link=logit)));
			},
			COX =
			{
			  fittedModel <- try(survival::coxph(model.formula,data=data,,na.action=na.exclude,model=TRUE,...));
			},
			{
	          fittedModel <- try(glm(model.formula,data=data,na.action=na.exclude,...));
#			  fittedModel <- try(speedglm::speedglm(model.formula,data=data,na.action=na.exclude,family=binomial(link=logit)));
			}
		)
		if (!inherits(fittedModel, "try-error"))
		{
			s <- is.na(fittedModel$coefficients);
			if (!any(s))
			{		
				if ((max(fittedModel$coefficients)>1.0e10) || (min(fittedModel$coefficients)< -1.0e10))
				{
					class(fittedModel) <- c(class(fittedModel),"try-error");
				}
			}
			else
			{
				class(fittedModel) <- c(class(fittedModel),"try-error");
			}
		}
	}
	else
	{	
		modelMat <- model.matrix(model.formula,data);
	    if (type=="COX")
		{
			response <- as.matrix(data[,all.vars(model.formula)[1:2]]);
		}
		else
		{
			M <- model.frame(model.formula, data)
			response <- as.matrix(cbind(M[,1],M[,1]));
		}
		fittedModel <-.Call("modelFittingCpp",response,modelMat,type);
		class(fittedModel) <- "fitFRESA";
		fittedModel$estimations <- as.vector(fittedModel$coefficients);
		fittedModel$formula <- model.formula;
		if (type=="COX")
		{
			fittedModel$coefficients <- fittedModel$estimations[1:(length(fittedModel$estimations)/2)];
		}
		names(fittedModel$coefficients) <- colnames(modelMat);
		
		fittedModel$type <- type;
		fittedModel$call <- match.call();
		fittedModel$terms <- terms(model.formula);
		s <- is.na(fittedModel$coefficients);
		if (!any(s))
		{		
			if ((max(fittedModel$coefficients)>1.0e10) || (min(fittedModel$coefficients)< -1.0e10))
			{
				class(fittedModel) <- c(class(fittedModel),"try-error");
			}
		}
		else
		{
			class(fittedModel) <- c(class(fittedModel),"try-error");
		}
	}
    return (fittedModel)
}
