ice = function(object, X, y,
		predictor, predictfcn, 
		verbose = TRUE, frac_to_build = 1, indices_to_build = NULL, 
		num_grid_pts, logodds = FALSE, probit = FALSE, ...){

	MAX_NUM_UNIQUE_PTS_NOMINAL = 5

	#check for factor
	if (class(X[, predictor]) == "factor" || class(X[, predictor]) == "character"){
		stop("ICE does not support factor attributes")
	}
	
	if(!is.numeric(frac_to_build) || frac_to_build > 1 || frac_to_build < 0 ){
		stop("frac_to_build must be in (0, 1]")
	}
  
  if(!missing(y) && class(y) == "factor"){
    stop("Do not pass y when it is categorical variable.")
  }
  
  if (logodds && probit){
	  stop("You must employ either logodds OR probit but not both.")
  }
  

	######## (1) check inputs
	# (a) check for valid prediction routine...
	if (!missing(predictfcn)){
		fcn_args = names(formals(predictfcn))
		if (!("object" %in% fcn_args && "newdata" %in% fcn_args)){
			stop("predictfcn must have 'object' and 'newdata' arguments")
		} else {
			use_generic = FALSE
		}
	} else { 
		#check for default prediction associated with class(object)
		#this is harder than it should be because some classes have class(object)
		#return a vector (randomForest, for instance).
		classvec = class(object) #may have multiple classes
		found_predict = FALSE
		i = 1
		for (i in 1 : length(classvec)){
			if (length(grep("predict", methods(class = classvec[i]))) > 0){
				found_predict = TRUE
				break
			}
		}
		if (!found_predict){
			stop("No generic predict method found for this object.")
		} else {
			use_generic = TRUE
		}
	}

	
	######## (2)
	N = nrow(X)
	# grid points
	#now create xj-to-predict values
	xj = X[, predictor]  #fix so predictor can be given as a name. 
	grid_pts = sort(X[, predictor])	

	# there are 3 cases: frac_to_build specificed, indices specified, or nothing specified.
	# 1: check fraction to build
	if (frac_to_build < 1){
		# we don't sample randomly -- we ensure uniform sampling across
		# quantiles of xj so as to not leave out portions of the dist'n of x.
		order_xj = order(xj)
		X = X[order_xj, , drop = FALSE]  #ordered by column xj 	
		nskip = round(1 / frac_to_build)
		xj_sorted_indices_to_build = seq(1, N, by = nskip)
		X = X[xj_sorted_indices_to_build, , drop = FALSE]
		xj = X[, predictor]
		grid_pts = sort(xj)
		indices_to_build = order_xj[xj_sorted_indices_to_build]
	} else { 
		
	  #2: indices specified:
	  if (!missing(indices_to_build)){
	  	#extract the indicies the user asks for first, THEN order the remaining
	  	#X matrix by the xj's left in the sub-sample defined by the user.
		if (frac_to_build < 1){
			stop("\"frac_to_build\" and \"indices_to_build\" cannot both be specified simultaneously")
		}
		X = X[indices_to_build, , drop = FALSE]
		xj = X[, predictor]
		order_xj = order(xj) 
	    #order the remaining by column xj
	    X = X[order_xj, ]  #ordered by column xj 	
		grid_pts = sort(xj)
		xj = X[, predictor] #end if for indices checking.
	  } else { #3: nothing specified, so just re-order by xj
		X = X[order(xj), , drop = FALSE]
		xj = X[, predictor]
	  }	
	}
	grid_pts = unique(grid_pts)
	num_unique_pts = length(grid_pts)
	
	
	#now handle less grid pts
	if (!missing(num_grid_pts)){
		if (num_grid_pts > num_unique_pts){
			warning(paste("the number of grid points specified,", num_grid_pts, "is larger than the number of unique values,", num_unique_pts, "(defaulting to number of unique values)"))	
		} else {
			grid_pts = grid_pts[round(seq(from = 1, to = num_unique_pts, length.out = num_grid_pts))]	
		}		
	}
	
	# generate partials
	if (use_generic){
		actual_prediction = predict(object, X, ...)
	} else {
		actual_prediction = predictfcn(object = object, newdata = X)
	}
	
	if (class(actual_prediction) == "factor"){
		stop("The predict function must return probabilities (not levels of a factor).")
	}
	if (logodds || probit){	
		min_pred = min(actual_prediction)
		max_pred = max(actual_prediction)
		
		#do some p_hat \in [0, 1] checking
		if (min_pred < 0){ 
			stop("the logodds option is on but predict returns values less than 0 (these should be probabilities!)")
		} else if (max_pred > 1){
			stop("the logodds option is on but predict returns values greater than 1 (these should be probabilities!)")
		}
		if (min_pred == 0){
			second_lowest = min(actual_prediction[actual_prediction > 0])
			if (is.na(second_lowest)){ 
				second_lowest = .0001
			}
			lowest_practical_prob = mean(c(second_lowest, 0))
			actual_prediction[actual_prediction == 0] = lowest_practical_prob
			#warning("At least one probability was predicted to be 0, ICEbox is using ", lowest_practical_prob, " for the value(s) instead.")
		}
		if (max_pred == 1){
			second_highest = max(actual_prediction[actual_prediction < 1])
			if (is.na(second_highest)){ 
				second_highest = .9999
			}
			highest_practical_prob = mean(c(second_highest, 1)) #arbitrarily, halfway between 1 and second_highest
			actual_prediction[actual_prediction == 1] = highest_practical_prob
			#warning("At least one probability was predicted to be 1, ICEbox is using ", highest_practical_prob, " for the value(s) instead.")
		}
		
		if (logodds){
			#centered logit formula
			actual_prediction = log(actual_prediction) - (1 / 2) * (log(actual_prediction) + log(1 - actual_prediction))
		} else if (probit){
			actual_prediction = qnorm(actual_prediction)
		}
		
	}
	
	ice_curves = matrix(NA, nrow = nrow(X), ncol = length(grid_pts))
	colnames(ice_curves) = grid_pts

	#Compute actual pdp. Note that this is averaged over the observations
	#we sample, so this might be different from the 'true' pdp if frac_to_build < 0.
	xvec_temp = X[ ,predictor]  #actual setting of this column after ordering.
	for (t in 1 : length(grid_pts)){
		X[, predictor] = grid_pts[t]
		if (use_generic){
			ice_curves[, t] = predict(object, X, ...)
		}
		else{
			ice_curves[, t] = predictfcn(object = object, newdata = X)
		}
		
		if(verbose){cat(".")}			
	}
	#return X to its original state.
	X[, predictor] = xvec_temp

	#do logit if necessary
	if (logodds || probit){
		#prevent log(0) error
		min_val = min(ice_curves)
		max_val = max(ice_curves)
		if (min_val < 0){
			stop("logodds is TRUE but predict returns negative values (these should be probabilities!)")
		} 
		if (min_val == 0){
			second_lowest = min(ice_curves[ice_curves > 0])
			if (is.na(second_lowest)){ 
				second_lowest = .0001 #arbitrary epsilon value
			}
			lowest_practical_prob = mean(c(second_lowest, 0)) #arbitrarily, halfway between 0 and second_lowest
			ice_curves[(ice_curves == 0)] = lowest_practical_prob
			warning("At least one probability was predicted to be 0, ICEbox is using ", lowest_practical_prob, " for the value(s) instead.")
		} 
		if (max_val == 1){
			second_highest = max(ice_curves[ice_curves < 1])
			if (is.na(second_highest)){ 
				second_highest = .9999 #arbitrary 1 - epsilon value
			} 
			highest_practical_prob = mean(c(second_highest, 1)) #arbitrarily, halfway between 1 and second_highest
			ice_curves[(ice_curves == 1)] = highest_practical_prob
			warning("At least one probability was predicted to be 1, ICEbox is using ", highest_practical_prob, " for the value(s) instead.")
		}
		if (logodds){
			#centered logit formula
			ice_curves = log(ice_curves) - (1 / 2) * (log(ice_curves) + log(1 - ice_curves)) 
		} else if (probit){
			ice_curves = qnorm(ice_curves)
		}		

	}
	if (verbose){cat("\n")}
	
	if (class(predictor) != "character"){
		xlab = paste("x", predictor, sep = "_")  #x_1, x_2 etc.
	} else {
		xlab = predictor #the actual name of the feature.
	}
	
	if (num_unique_pts <= MAX_NUM_UNIQUE_PTS_NOMINAL){
		nominal_axis = TRUE
	} else {
		nominal_axis = FALSE
	}	

	range_y = NULL
	sd_y = NULL
	if(!missing(y)){
		range_y = max(y) - min(y)
		sd_y = sd(y)
	} else if (!logodds && !probit){
		range_y = (max(ice_curves) - min(ice_curves))
		sd_y = sd(actual_prediction)
		cat("y not passed, so range_y is range of ice curves and sd_y is sd of predictions on real observations\n")
	}

	#compute friedman's pdp:
	pdp = apply(ice_curves, 2, mean) # pdp = average over the columns
  	if(missing(predictfcn)){
	    predictfcn=NULL
	}
	
	ice_obj = list(ice_curves = ice_curves, gridpts = grid_pts, predictor = predictor, xj = xj, actual_prediction = actual_prediction, 
			logodds = logodds, probit = probit, xlab = xlab, nominal_axis = nominal_axis, range_y = range_y, sd_y = sd_y, Xice = X, pdp = pdp,
			indices_to_build = indices_to_build, frac_to_build = frac_to_build, predictfcn = predictfcn) 
	class(ice_obj) = "ice"
		
	invisible(ice_obj)
}
