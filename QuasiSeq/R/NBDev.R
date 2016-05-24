NBDev <- function(counts, design, log.offset, nb.disp, print.progress = TRUE, bias.fold.tolerance=1.10) {

    n <- ncol(counts)
    
    if (is.null(log.offset)) 
        log.offset <- rep(0, ncol(counts))
    est.offset <- exp(log.offset)
    
    deviance.vector <- rep(NA, nrow(counts))
    means <- matrix(NA, nrow(counts), ncol(counts))
    parms <- matrix(NA, nrow(counts), ncol(design))
	design.df=as.data.frame(design)
	glm.ctrl=glm.control(epsilon = 1e-08, maxit = 1500L, trace = FALSE)
	fbrNBglm.ctrl=fbrNBglm.control(coefOnly=TRUE, infoParms=list(j=1,k=1,m=1), maxit=1500L, tol=1e-8, standardizeX=TRUE)
	nbLogFamily=negbin("log", 1)
	logFcCorrection=abs(log(bias.fold.tolerance))
	log10=log(10)
	
    ### For each gene and given design matrix, fit glm using provided negative binomial dispersion estimate
    for (gn in 1:nrow(counts)) {
        ### If wanted, provide running progress update (eventually once every 5000 genes)
        if (gn %in% c(2, 10, 100, 500, 1000, 2500, 5000 * (1:200)) & print.progress) 
            print(paste("Analyzing Gene #", gn))
        
        #### For 2000 Fly genes, glm takes roughly 9 seconds (optim took roughly 21 seconds)
		glm.fitted <- glmsolve(formula = counts[gn, ] ~ . - 1 + offset(log.offset), family = update.fbrNBfamily(nbLogFamily, overdisp=nb.disp[gn]), data = design.df, control = glm.ctrl, x=TRUE) 
			
        parms[gn, ] <- glm.fitted$coefficients
        
        ### Save deviance (used to compute LRT for comparing models and also deviance dispersion estimator)
        deviance.vector[gn] <- glm.fitted$deviance
        
        ### When a group has only zero counts, the MLE doesn't exist.  For these genes, we use the method Kosmidis & Firth(2009)
        ### to moderate the parameter estimates.
        ### For 2000 Fly genes, fbr takes roughly 41 seconds
        #if (any(counts[gn, ] == 0) && any(abs(coefficients(glm.fitted)) > 3)) {
		tmp.bias = coef(glm.fitted, type='bias')
		tmp.nonNA=which(!is.na(tmp.bias) & !is.na(glm.fitted$coefficients))
		this.x=model.matrix(glm.fitted)
		if(any(abs(this.x[,tmp.nonNA,drop=FALSE]%*%tmp.bias[tmp.nonNA]) > logFcCorrection)) { 
			fbrNBglm.ctrl$start = glm.fitted$coefficients - pmax.int(-log10, pmin.int(log10, tmp.bias))
			wt.idx=which(glm.fitted$weights>0)
            fbrglm.fit <- fbrNBglm.fit(x=this.x[wt.idx,,drop=FALSE], y=counts[gn, wt.idx],offset=log.offset[wt.idx], family=nbLogFamily, control = fbrNBglm.ctrl)
			## note: As fbrNBglm.fit is only called right after a update.fbrfamily call, the overdisp parameter is already up-to-date
            parms[gn, ] <- fbrglm.fit
			tmp.nonNA = which(!is.na(fbrglm.fit))
        }
		
        ### Save optimized means (used in Pearson's dispersion estimator)
		### NOTE: This line was originally before the Firth bias correction block. 
		means[gn, ] <- as.vector(exp(design[,tmp.nonNA,drop=FALSE] %*% parms[gn, tmp.nonNA]) * est.offset)
    }
    return(list(dev = deviance.vector, means = means, parms = parms))
} 



