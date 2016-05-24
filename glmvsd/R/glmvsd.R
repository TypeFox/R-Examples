glmvsd <- function(x, y, n_train = ceiling(n/2), no_rep = 100, 
					n_train_bound = n_train-2, n_bound = n-2,
					model_check, psi = 1, 
					family = c("gaussian", "binomial"), method = c("union", "customize"),
                   	candidate_models, weight_type = c("BIC", "AIC", "ARM"), prior = TRUE,
					reduce_bias = FALSE) {
    # check data and parameter
    family <- match.arg(family)
    method <- match.arg(method)
    weight_type <- match.arg(weight_type)
    y <- drop(y)
    y <- as.numeric(y)
    x <- as.matrix(x)
    p <- NCOL(x)
    n <- length(y)
    if (family == "binomial") {
      if (!all(y %in% c(0, 1))) 
        stop("There can only be 0 or 1 in y when using binomial family")
    }
    if (n != NROW(x)) 
        stop("x and y have different number of observations")
    if (n_train >= n) 
        stop("Training size must be less than the number of observations")
    if (missing(model_check)) 
        stop("User must provide a base model.")
	if(is.vector(model_check)) model_check <- matrix(model_check,nrow=1)
    # use union option to compute candidate models
    if (method == "union") {
      if (family == "gaussian") 
        candidate_models <- gaussianfit(x, y)
      if (family == "binomial") 
        candidate_models <- binomialfit(x, y)
    }
    if (method == "customize") {
      if (missing(candidate_models)) 
        stop("Users must supply a candidate model.")
      if (!is.matrix(candidate_models)) 
        stop("Supplied model must be a matrix.")
      if (NCOL(candidate_models) != NCOL(x)) 
        stop("Number of variables in candidate model and x does not match.")
      if (!all(as.numeric(candidate_models) %in% c(0, 1))) 
        stop("There can only be 0 or 1 in candidate_models")
    }
    # clean the candidate models mk
    candidate_models <- unique(candidate_models)
    rownames(candidate_models) <- NULL
    candidate_models <- candidate_models[order(rowSums(candidate_models)), ]
    if (weight_type == "ARM") {
      candidate_models <- candidate_models[rowSums(candidate_models) < n_train_bound, ]
    }else{
	  candidate_models <- candidate_models[rowSums(candidate_models) < n_bound, ]
    }
    # compute weights
    if (family == "gaussian") {
      if (weight_type == "ARM") {
        fit <- lsARM(x = x, y = y, candidate_models = candidate_models, 
                     n_train = n_train, no_rep = no_rep, psi = psi, prior = prior)
      }
      if (weight_type == "AIC" | weight_type == "BIC") {
        fit <- lsIC(x = x, y = y, candidate_models = candidate_models, 
                    psi = psi, prior = prior, type = weight_type)
      }
    }
    if (family == "binomial") {
      if (weight_type == "ARM") {
        fit <- logitARM(x = x, y = y, candidate_models = candidate_models, 
                        n_train = n_train, no_rep = no_rep, psi = psi, 
						prior = prior, reduce_bias = reduce_bias)
      }
      if (weight_type == "AIC" | weight_type == "BIC") {
        fit <- logitIC(x = x, y = y, candidate_models = candidate_models, 
						psi = psi, prior = prior, type = weight_type, 
						reduce_bias = reduce_bias)
      }
    }
    weight <- fit$weight
	  # compute VSD etc
	  # initialization
	  no_checkmod <- NROW(model_check)
	  VSD <- rep(NA, no_checkmod)
	  VSD_minus <- rep(NA, no_checkmod)
	  VSD_plus <- rep(NA, no_checkmod)
	  Fmeasure <- rep(NA, no_checkmod)
	  Gmeasure <- rep(NA, no_checkmod)
	  # size of m0
	  model_check_size <- rowSums(model_check)
	  # size of mk
	  candidate_models_size <- rowSums(candidate_models)
	  # start the loop
	  for (mindex in seq(no_checkmod)) {
	    # compare m0 and mk
	    TMP_matrix <- sweep(candidate_models, MARGIN = 2, model_check[mindex, ], "-")
	    diff <- rowSums(abs(TMP_matrix))
	    diff_plus <- rowSums(TMP_matrix == 1)
	    diff_minus <- rowSums(TMP_matrix == -1)
	    # compute VSD and VSD plus and minus
	    VSD[mindex] <- sum(weight*diff)  # glmvsd value
	    VSD_plus[mindex] <- sum(weight*diff_plus)  # false negative 
	    VSD_minus[mindex] <- sum(weight*diff_minus)  # false positive
	    # compute F measure and G measure using precision and recall
	    Prcision <- (model_check_size[mindex]-diff_minus)/model_check_size[mindex]
	    Recall <- (model_check_size[mindex]-diff_minus)/candidate_models_size
	    Fmeasure_tmp <- 2*(Prcision*Recall)/(Prcision+Recall)
	    Gmeasure_tmp <- sqrt(Prcision*Recall)
	    Fmeasure[mindex] <- sum(weight*Fmeasure_tmp, na.rm = TRUE)
	    Gmeasure[mindex] <- sum(weight*Gmeasure_tmp, na.rm = TRUE)
	  }
    # output 
    object <- list(candidate_models_cleaned = candidate_models, VSD = VSD, VSD_minus = VSD_minus, VSD_plus = VSD_plus,  
	    			Fmeasure = Fmeasure, Gmeasure = Gmeasure,
                   weight = weight)
    class(object) <- "glmvsd"
    object
}