run_model_on_left_out_record_results_and_cleanup = function(Xy, 
		leave_outs_to_be_predicted, 
		train_on_all_except_these, 
		model_string, 
		predict_string, 
		cleanup_mod_function, 
		y_higher_is_better,
		full_verbose = FALSE, 
		verbose = FALSE,
		...){
	
	#the left one out matrix has n-1 rows and will be considered the "training data"
	Xyleft = Xy[-train_on_all_except_these, ]
	
	#build the model via the user-specified string
	mod = eval(parse(text = model_string)) #this function makes use of the "Xyleft" object
	
	#pull out the record of the left-one-out subject
	obs_left_out = Xy[leave_outs_to_be_predicted, 1 : (ncol(Xy) - 1)]
	
	#also take note of what actually happened to this subject in the experiment
	real_ys = Xy[leave_outs_to_be_predicted, ncol(Xy)]
	orig_trts = obs_left_out$treatment
	
	#now evaluate the left-one-out subject on the model for both his true treatment and his counterfactual
	obs_left_out$treatment = 0
	yhatTx0s = eval(parse(text = predict_string))
	obs_left_out$treatment = 1
	yhatTx1s = eval(parse(text = predict_string))
	
	#give the user some indication of progress if they want to see it
	if (full_verbose){
		cat("model #", leave_outs_to_be_predicted, "/", nrow(Xy), " yhatTx0:1 = ", round(yhatTx0s, 2), " : ", round(yhatTx1s, 2), "\n", sep = "")
	} else if (verbose){
		cat(".")
	}
	
	#if the models need to be cleaned up in some way, do it now before the next iteration of the leave-one-out
	if (!is.na(cleanup_mod_function)){
		eval(parse(text = cleanup_mod_function))
	}
	
	#tabulate the result for the prediction on this left one out model
	tabulate_results_for_left_one_out_subject(orig_trts, yhatTx0s, yhatTx1s, real_ys, y_higher_is_better)
}