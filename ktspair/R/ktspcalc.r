ktspcalc <- function (dat, grp, k = NULL, cross = 5, performance = FALSE, healthy = NULL, seed = NULL, display = TRUE, length = 40, med = FALSE) {

## Compute the k-TSP by using the function kts.pair(). If no parameter is inserted for k, a crossvalidation is performed.

	accuracy <- NULL
	accuracy_k <- NULL
	sensitivity <- NULL
	specificity <- NULL

	if(!is.null(seed)){set.seed(seed)}

	if(!is.null(healthy)){
		if(!(healthy %in% grp)) stop("The group called ",healthy," is not present in the variable grp \n")
		else performance <- TRUE
	}

	if(!is.null(k)){
		if(k>length){stop("The value of length has to be at least as big as k")}
		if(k%%2==0){stop("k needs to be an odd number")}
	}

	if(is.null(k) && cross > ncol(dat)){cat("The number of fold is too high, it should be between 1 and ", ncol(dat)," \n")}

	if (class(dat) == "ExpressionSet"){
        genenames <- as.character(1:dim(exprs(dat))[1])
        if (!is.null(featureNames(dat))) {
			genenames <- featureNames(dat)
        }
        if (length(grp) == 1) {
			labels <- as.character(unique(pData(dat)[, grp])[order(unique(grp))])
			if(!is.null(healthy)){healthy.position <- which((pData(dat)[, grp])==healthy)[1]}
			grp <- make.consecutive.int(pData(dat)[, grp])## Transform the variable grp into a binary vector.
			if(!is.null(healthy)){healthy <- grp[healthy.position]}
			dat <- exprs(dat)
			rownames(dat) <- genenames
				if (max(grp) != 1) {
					stop("TSPs can only be calculated for variables with two classes")
				}
			if(is.null(k)){## If k was not inserted by the user.
				cv <- cv(dat, grp, cross = cross, display = display, length = length, med = med, healthy = healthy)
				k <- cv$k
				if(!is.null(length)){length <- 3*k}
				accuracy_k <- cv$accuracy_k
				accuracy <- cv$accuracy
				if(!is.null(healthy)){
					sensitivity <- cv$sensitivity
					specificity <- cv$specificity
				}
			}
			if(performance && !is.null(k)){## If k was not inserted by the user and the performances of the model are required.
				cv <- cv(dat, grp, cross = cross, display = display, length = length, med = med, healthy = healthy)
				accuracy_k <- cv$accuracy_k
				accuracy <- cv$accuracy
				if(!is.null(healthy)){
					sensitivity <- cv$sensitivity
					specificity <- cv$specificity
				}
			}
			ktsp <- kts.pair(dat, grp, k, display = display, length = length, med = med)
			ktsp$labels <- labels
			ktsp$accuracy_k <- accuracy_k
			ktsp$accuracy <- accuracy
			if(!is.null(healthy)){
				ktsp$sensitivity <- sensitivity
				ktsp$specificity <- specificity
			}
			return(ktsp)
		}
    
		else{
			labels <- as.character(unique(grp)[order(unique(grp))])
			if(!is.null(healthy)){healthy.position <- which(grp==healthy)[1]}
			grp <- make.consecutive.int(grp)
			if(!is.null(healthy)){healthy <- grp[healthy.position]}
			dat <- exprs(dat)
			rownames(dat) <- genenames
			if (max(grp) != 1) {
				stop("TSPs can only be calculated for variables with two classes")
			}
			if(!is.null(healthy)){healthy <- which(labels==healthy)-1}
			if(is.null(k)){
				cv <- cv(dat, grp, cross = cross, display = display, length = length, med = med, healthy = healthy)
				k <- cv$k
				if(!is.null(length)){length <- 3*k}
				accuracy_k <- cv$accuracy_k
				accuracy <- cv$accuracy
				if(!is.null(healthy)){
					sensitivity <- cv$sensitivity
					specificity <- cv$specificity
				}
			}
			if(performance && !is.null(k)){
				cv <- cv(dat, grp, cross = cross, display = display, length = length, med = med, healthy = healthy)
				accuracy_k <- cv$accuracy_k
				accuracy <- cv$accuracy
				if(!is.null(healthy)){
					sensitivity <- cv$sensitivity
					specificity <- cv$specificity
				}
			}
			ktsp <- kts.pair(dat, grp, k, display = display, length = length, med = med)
			ktsp$labels <- labels
			ktsp$accuracy_k <- accuracy_k
			ktsp$accuracy <- accuracy
			if(!is.null(healthy)){
				ktsp$sensitivity <- sensitivity
				ktsp$specificity <- specificity
			}
			return(ktsp)
		}
	}
	else {
        labels <- as.character(unique(grp)[order(unique(grp))])
	if(!is.null(healthy)){healthy.position <- which(grp==healthy)[1]}
	grp <- make.consecutive.int(grp)
	if(!is.null(healthy)){healthy <- grp[healthy.position]}
        if (is.null(rownames(dat))) {
			rownames(dat) <- as.character(1:dim(dat)[1])
        }
        if (max(grp) != 1) {
			stop("TSPs can only be calculated for variables with two classes")
        }		
		if(!is.null(healthy)){healthy <- which(labels==healthy)-1}
		if(is.null(k)){
			cv <- cv(dat, grp, cross = cross, display = display, length = length, med = med, healthy = healthy)
			k <- cv$k
			accuracy_k <- cv$accuracy_k
			accuracy <- cv$accuracy
			if(!is.null(healthy)){
				sensitivity <- cv$sensitivity
				specificity <- cv$specificity
			}
		}
		if(performance && !is.null(k)){
			cv <- cv(dat, grp, cross = cross, display = display, length = length, med = med, healthy = healthy)
			accuracy_k <- cv$accuracy_k
			accuracy <- cv$accuracy
			if(!is.null(healthy)){
				sensitivity <- cv$sensitivity
				specificity <- cv$specificity
			}
		}
        ktsp <- kts.pair(dat, grp, k, display = display, length = length, med = med)
        ktsp$labels <- labels
	ktsp$accuracy_k <- accuracy_k
	ktsp$accuracy <- accuracy
	if(!is.null(healthy)){
		ktsp$sensitivity <- sensitivity
		ktsp$specificity <- specificity
	}
        return(ktsp) 
	}
}
