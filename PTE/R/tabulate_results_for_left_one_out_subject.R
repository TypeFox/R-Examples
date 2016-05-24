tabulate_results_for_left_one_out_subject = function(orig_trts, yhatTx0s, yhatTx1s, real_ys, y_higher_is_better){
	res = matrix(NA, nrow = length(orig_trts), ncol = 5)
	
	for (i_left_out in 1 : length(orig_trts)){
		est_true = ifelse(orig_trts[i_left_out] == 0, yhatTx0s[i_left_out], yhatTx1s[i_left_out])
		est_counterfactual = ifelse(orig_trts[i_left_out] == 0, yhatTx1s[i_left_out], yhatTx0s[i_left_out])
		if (y_higher_is_better){
			optimal = as.numeric(est_true > est_counterfactual)
		} else {
			optimal = as.numeric(est_true < est_counterfactual)
		}
		res[i_left_out, ] = c(est_true, est_counterfactual, as.numeric(orig_trts[i_left_out] == 0), optimal, real_ys[i_left_out])		
	}	
	res
}