#' @method summary bootstrapValidation_Bin

summary.bootstrapValidation_Bin <- 
function(object, ...) 
{


	 digits = 3 
	conf.int = 0.95
	cilow = (1.0-conf.int)/2;
	cihig = 1.0-cilow;

	ciAccuracy <- as.vector(quantile(object$train.accuracy, probs = c(cilow,0.5, cihig), na.rm = TRUE,names = FALSE, type = 7));
	ciSensitivity <- as.vector(quantile(object$train.sensitivity, probs = c(cilow,0.5, cihig), na.rm = TRUE,names = FALSE, type = 7));
	ciSpecificity <- as.vector(quantile(object$train.specificity, probs = c(cilow,0.5, cihig), na.rm = TRUE,names = FALSE, type = 7));

	citROC <- as.vector(quantile(object$train.ROCAUC, probs = c(cilow,0.5, cihig), na.rm = TRUE,names = FALSE, type = 7));
	
		
	bootAuc <- pROC::roc( object$outcome, object$boot.model$linear.predictors,plot=FALSE,ci=TRUE,boot.n=object$loops);
	
    cat("\nModel Cross-Validation with Improvement in Predicted Probability\n\n")
    cat("Number of Cases:", sum(object$outcome), "\t Number of Controls", sum(object$outcome==0), "\n\n")
    cat("Number of Bootstraps:", length(object$train.ROCAUC), "\t Sampled Fraction", object$fraction, "\n\n")

	performance <- vector();
	
	cat(l1 <- sprintf("Blind    Accuracy: %8.3f : Bootstrapped    Accuracy: %8.3f (%8.3f to %8.3f) \n",object$blind.accuracy,ciAccuracy[2],ciAccuracy[1], ciAccuracy[3]));
	performance <- append(performance,l1);
	cat(l1 <- sprintf("Blind Sensitivity: %8.3f : Bootstrapped Sensitivity: %8.3f (%8.3f to %8.3f) \n",object$blind.sensitivity,ciSensitivity[2],ciSensitivity[1], ciSensitivity[3]));
	performance <- append(performance,l1);
	cat(l1 <- sprintf("Blind Specificity: %8.3f : Bootstrapped Specificity: %8.3f (%8.3f to %8.3f) \n",object$blind.specificity,ciSpecificity[2],ciSpecificity[1], ciSpecificity[3]));
	performance <- append(performance,l1);
	cat(l1 <- sprintf("Blind     ROC AUC: %8.3f : Bootstrapped     ROC AUC: %8.3f (%8.3f to %8.3f) \n",object$blind.ROCAUC$auc,object$boot.ROCAUC$auc,bootAuc$ci[1],bootAuc$ci[3]));
	performance <- append(performance,l1);
	cat(l1 <- sprintf("Blind     ROC AUC: %8.3f : ModelBootstrap   ROC AUC: %8.3f (%8.3f to %8.3f) \n\n",object$blind.ROCAUC$auc,citROC[2],citROC[1],citROC[3]));
	performance <- append(performance,l1);
	
	performance.table <- rbind(c(object$blind.accuracy,ciAccuracy[2],ciAccuracy[1], ciAccuracy[3]));
	performance.table <- rbind(performance.table,c(object$blind.sensitivity,ciSensitivity[2],ciSensitivity[1], ciSensitivity[3]));
	performance.table <- rbind(performance.table,c(object$blind.specificity,ciSpecificity[2],ciSpecificity[1], ciSpecificity[3]));
	performance.table <- rbind(performance.table,c(object$blind.ROCAUC$auc,object$boot.ROCAUC$auc,bootAuc$ci[1],bootAuc$ci[3]));
	performance.table <- rbind(performance.table,c(object$blind.ROCAUC$auc,citROC[2],citROC[1],citROC[3]));
	colnames(performance.table) <- c("Blind","Train","LCI","UCI")
	rownames(performance.table) <- c("Accuracy","Sensitivity","Specificity","ROCAUC 1","ROCAUC 2")
	
	
	smry <- summary(object$boot.model, ...);
	meancoef <- colMeans(object$s.coef,na.rm = TRUE);
	lowci <- vector();
	topci <- vector();
	for ( i in 1:length(meancoef))
	{
		ci <- as.vector(quantile(object$s.coef[,i], probs = c(cilow, cihig), na.rm = TRUE,names = FALSE, type = 7));
		lowci <- append(lowci,ci[1]);
		topci <- append(topci,ci[2]);
	}
	zidimedian <- vector();
	idimedian <- vector();
	idilowci <- vector();
	iditopci <- vector();
	znrimedian <- vector();
	nrimedian <- vector();
	nrilowci <- vector();
	nritopci <- vector();

	classlen=length(class(object$boot.model))
	if (substr(class(object$boot.model)[classlen], 1, 2) == "co")
	{
		startlist = 1;
	}
	else
	{
		idimedian <- append(idimedian,0);
		zidimedian <- append(zidimedian,0);
		idilowci <- append(idilowci,0);
		iditopci <- append(iditopci,0);
		nrimedian <- append(nrimedian,0);
		znrimedian <- append(znrimedian,0);
		nrilowci <- append(nrilowci,0);
		nritopci <- append(nritopci,0);
		startlist = 2;
	}

	
	for ( i in startlist:length(meancoef))
	{
		j=i+1-startlist;
		ci <- as.vector(quantile(object$IDIs[,j], probs = c(cilow, 0.5, cihig), na.rm = TRUE,names = FALSE, type = 7));
		idimedian <- append(idimedian,ci[2]);
		idilowci <- append(idilowci,ci[1]);
		iditopci <- append(iditopci,ci[3]);		
		ci <- as.vector(quantile(object$z.IDIs[,j], probs = c(cilow, 0.5, cihig), na.rm = TRUE,names = FALSE, type = 7));
		zidimedian <- append(zidimedian,ci[2]);

		ci <- as.vector(quantile(object$NRIs[,j], probs = c(cilow, 0.5, cihig), na.rm = TRUE,names = FALSE, type = 7));
		nrimedian <- append(nrimedian,ci[2]);
		nrilowci <- append(nrilowci,ci[1]);
		nritopci <- append(nritopci,ci[3]);		
		ci <- as.vector(quantile(object$z.NRIs[,j], probs = c(cilow, 0.5, cihig), na.rm = TRUE,names = FALSE, type = 7));
		znrimedian <- append(znrimedian,ci[2]);

	}


	p <- smry$coefficients;
	p <- cbind(p,lowci);
	p <- cbind(p,topci);
	p <- cbind(p,idimedian);
	p <- cbind(p,idilowci);
	p <- cbind(p,iditopci);
	p <- cbind(p,zidimedian);

	p <- cbind(p,nrimedian);
	p <- cbind(p,nrilowci);
	p <- cbind(p,nritopci);
	p <- cbind(p,znrimedian);

	cnames <- colnames(smry$coefficients);
	cnames <- append(cnames,c("Low CI","High CI","Median IDI","Low IDI","High IDI","z IDI","Median NRI","Low NRI","High NRI","z NRI"));
	colnames (p) <- cnames;
	rownames(p) <- colnames(object$s.coef);
	print(p, digits = digits);
	
	result <- list(performance = performance,
	summary = smry,
	coef = p, 
	performance.table = performance.table);
	
	return (result);
        
}
