# Common methods for different model classes

# Method analogous to model.frame.lm
getModelFrame.default <- function(model){model.frame(model)}
getModelFrame.lme <- function(model){
	mt <- model$terms
	mf <- model$data[names(attr(mt,"dataClasses"))]
	attr(mf,"terms") <- mt
	return(mf)
}
getModelFrame <- function(model){UseMethod("getModelFrame")}

# Method analogous to model.matrix.lm
getModelMatrix.default <- function(model){model.matrix(model)}
getModelMatrix.lme <- function(model){model.matrix(terms(model),getModelFrame(model))}
getModelMatrix <- function(model){UseMethod("getModelMatrix")}

# Method analogous to coef.lm
getCoef.default <- function(model){coef(model)}
getCoef.merMod <- getCoef.mer <- getCoef.lme <- function(model){lme4::fixef(model)}
getCoef.lme <- function(model){nlme::fixef(model)}
getCoef <- function(model){UseMethod("getCoef")}

# Method for obtaning factor levels from model
getXLevels.lm <- function(model){model$xlevels}
getXLevels.lme <- function(model){lapply(model$contrasts, rownames)}
# Extract factors from a list (use it in the default method)
getXLevels.list <- function(fr){
	are.factors <- sapply(fr, is.factor)
	lapply(fr[are.factors], "levels")
}
getXLevels.default <- function(model){
	predictors <- getPredictors(model)
	mf <- model.frame(model)[predictors]
	getXLevels.list(mf)
}
getXLevels <- function(model){UseMethod("getXLevels")}

# Method for obtaining contrasts
getContrasts <- function(model){attr(getModelMatrix(model),"contrasts")}

# Method for defining the error family (in glm or glmm)
getFamily.default <- function(model){NULL}
getFamily.merMod <- getFamily.glm <- function(model){family(model)}
getFamily.mer <- function(model){
	if (length(model@muEta > 0)){ # isGLMM
		# Get family from function call
		fam <- getCall(model)$family
		# If it was a character string, try get the appropriate function
		if (is.character(fam)){
			tryCatch(
				fam <- get(fam, mode = "function"),
				error=function(e) fam <- NULL
			)
		}
		return(eval(fam)())
	}else return(NULL)
}
getFamily <- function(model){UseMethod("getFamily")}

# Method for obtaining a character vector of the model predictors,
# excluding offsets
getPredictors <- function(model){
	modterms <- terms(model)
	model.variables <- as.character(attr(modterms,"variables"))[-1]
	model.variables[-c(attr(modterms,"response"), attr(modterms,"offset"))]
}

# get offset variable, if any
getOffset.lm <- getOffset.glm <- function(model){model$offset}
getOffset.mer <- function(model){model@offset}
getOffset.merMod <- function(model){lme4::getME(model,"offset")}
getOffset.default <- function(model){
	modterms <- terms(model)
	if (!is.null(offpos <- attr(modterms, "offset"))){
		with(eval(getCall(model)$data), eval(languageEl(attr(terms(model), "variables"), offpos+1L)))
	}else NULL
}
getOffset <- function(model){UseMethod("getOffset")}
