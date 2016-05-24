create_raw_results_matrix = function(n){
	raw_results = as.data.frame(matrix(NA, nrow = n, ncol = 5))
	colnames(raw_results) = c("est_true", "est_counterfactual", "trt0", "optimal", "real_y")
	raw_results
}