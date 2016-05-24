


cv <-function (dat, grp, cross = 5, display = FALSE, length = 40, seed = NULL, med = FALSE, healthy = NULL){## Perform, by defalut, a 5 fold crossvalidation on the dataset.
	
	if(!is.null(seed)){set.seed(seed)}

	if(!is.null(healthy)){if(!(healthy %in% grp)) stop("The group called ",healthy," is not present in the variable grp \n")}
	if(is.null(healthy)){healthy <- grp[1]}
	n <- ncol(dat)
	m <- nrow(dat)
	bound <- m/4
	accuracy <- numeric(length=5)
	sensitivity <- numeric(length=5)
	specificity <- numeric(length=5)
	set <- sample(rep(as.vector(1:cross),ceiling(dim(dat)[2]/cross)))[1:n]## Used to seperate the dataset into random partitions

	for (j in 0:min(4,bound)){
		correct_pred <- numeric(length=cross)
		correct_sensitivity <- numeric(length=cross)
		correct_specificity <- numeric(length=cross) 
		possible_fold<-0
		number <-c()
		for (i in 1:cross){
			ktsp <- ktspcalc(dat[,set!=i], grp[set!=i], k = 2 * j + 1, display = display, length = length, med = med)
			if(ktsp$k == (2*j+1)){
				dat_cross <- dat[,set==i]
				grp_cross <- grp[set==i]
				sensitivity_pred <- predict(ktsp, dat_cross[,grp_cross!=healthy], display=display) 
				specificity_pred <- predict(ktsp, dat_cross[,grp_cross==healthy], display=display)
				sensitivity_possible <- length(which(sensitivity_pred!=""))
				specificity_possible <- length(which(specificity_pred!=""))
				correct_sensitivity[i] <- sum(sensitivity_pred== grp_cross[grp_cross!=healthy])/sensitivity_possible
				correct_specificity[i] <- sum(specificity_pred== grp_cross[grp_cross==healthy])/specificity_possible
				prediction <- predict(ktsp, dat_cross, display=display) 
				pred_possible <- length(which(prediction!=""))
				correct_pred[i] <- sum(prediction== grp_cross)/pred_possible
				possible_fold <- possible_fold + 1
			}
		}
 	if(possible_fold>0){
		accuracy[j+1] <- mean(correct_pred[is.na(correct_pred)==FALSE])
		sensitivity[j+1] <- mean(correct_sensitivity[is.na(correct_sensitivity)==FALSE])
		specificity[j+1] <- mean(correct_specificity[is.na(correct_specificity)==FALSE])}## Compute the accuracy reached with different values of k
	if(display){
		if(possible_fold<0.5*cross){cat("Crossvalidation with k = ", 2*j+1 ," worked only for ", possible_fold," partitions \n")}
	}
	}
	cv <- list(k = 2 * which.max(accuracy)-1, accuracy_k = accuracy[which.max(accuracy)], accuracy = accuracy, sensitivity = sensitivity, specificity = specificity)
	class(cv) <- "cv"
	return(cv)	
}


