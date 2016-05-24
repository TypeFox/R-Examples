


cv2 <-function (indice, dat, grp, cross = 5, display = FALSE, length = 40, seed = NULL, med = FALSE, healthy = NULL){
	
	if(!is.null(seed)){set.seed(seed)}

	if(!is.null(healthy)){ if(!healthy %in% grp) stop("The group called ",healthy," is not present in the variable grp \n")}
	if(is.null(healthy)){healthy <- grp[1]}

	n <- ncol(dat)
	m <- nrow(dat)
	bound <- m/4-1/2

	accuracy <- NULL
	sensitivity <- NULL
	specificity <- NULL
	set <- sample(rep(as.vector(1:cross),ceiling(dim(dat)[2]/cross)))[1:n]## Used to seperate the dataset into random partitions


	correct_pred <- numeric(length=cross)
	correct_sensitivity <- numeric(length=cross)
	correct_specificity <- numeric(length=cross)

	number <-c()
	for (i in 1:cross){
		ktsp <- ktspcalc2(indice,dat[,set!=i], grp[set!=i],performance=FALSE, display = display, length = length, med = med)

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
	}


	accuracy <- mean(correct_pred[is.na(correct_pred)==FALSE])
	sensitivity <- mean(correct_sensitivity[is.na(correct_sensitivity)==FALSE])
	specificity <- mean(correct_specificity[is.na(correct_specificity)==FALSE])

	cv <- list(k = length(indice)/2, accuracy = accuracy, sensitivity = sensitivity, specificity = specificity)
	class(cv) <- "cv"
	return(cv)	
}


