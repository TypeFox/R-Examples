create_PTE_results_object = function(raw_results, y_higher_is_better){
	n = nrow(raw_results)
	
	#get real y's
	real_ys_for_summary = list()
	real_ys_for_summary[["received Trt 0, optimal Trt 0"]] = raw_results[raw_results$trt0 == 1 & raw_results$optimal == 1, "real_y"] #optimal
	real_ys_for_summary[["received Trt 0, optimal Trt 1"]] = raw_results[raw_results$trt0 == 1 & raw_results$optimal == 0, "real_y"] #not optimal
	real_ys_for_summary[["received Trt 1, optimal Trt 0"]] = raw_results[raw_results$trt0 == 0 & raw_results$optimal == 0, "real_y"] #not optimal
	real_ys_for_summary[["received Trt 1, optimal Trt 1"]] = raw_results[raw_results$trt0 == 0 & raw_results$optimal == 1, "real_y"] #optimal
	
	Ra = real_ys_for_summary[["received Trt 0, optimal Trt 0"]]
	Rb = real_ys_for_summary[["received Trt 0, optimal Trt 1"]]
	Rc = real_ys_for_summary[["received Trt 1, optimal Trt 0"]]
	Rd = real_ys_for_summary[["received Trt 1, optimal Trt 1"]]
	
	#four groups
	summary_square = matrix(NA, 3, 3)
	rownames(summary_square) = c("received Trt 0", "received Trt 1", "Avg Opt")
	colnames(summary_square) = c("optimal Trt 0", "optimal Trt 1", "Avg Trt")
	summary_square[1, 1] = mean(Ra, na.rm = TRUE)
	summary_square[1, 2] = mean(Rb, na.rm = TRUE)
	summary_square[2, 1] = mean(Rc, na.rm = TRUE)
	summary_square[2, 2] = mean(Rd, na.rm = TRUE)
	summary_square[1, 3] = mean(c(Ra, Rb), na.rm = TRUE)
	summary_square[2, 3] = mean(c(Rc, Rd), na.rm = TRUE)
	summary_square[3, 1] = mean(c(Ra, Rc), na.rm = TRUE)
	summary_square[3, 2] = mean(c(Rb, Rd), na.rm = TRUE)
	summary_square[3, 3] = mean(c(Ra, Rb, Rc, Rd), na.rm = TRUE)
	
	#four groups
	ns = matrix(NA, 3, 3)
	rownames(ns) = c("received Trt 0", "received Trt 1", "Totals")
	colnames(ns) = c("optimal Trt 0", "optimal Trt 1", "Totals")
	ns[1, 1] = length(Ra[!is.na(Ra)])
	ns[1, 2] = length(Rb[!is.na(Rb)])
	ns[2, 1] = length(Rc[!is.na(Rc)])
	ns[2, 2] = length(Rd[!is.na(Rd)])
	#now add em up to get tots
	ns[1, 3] = ns[1, 1] + ns[1, 2]
	ns[2, 3] = ns[2, 1] + ns[2, 2]
	ns[3, 1] = ns[1, 1] + ns[2, 1]
	ns[3, 2] = ns[1, 2] + ns[2, 2]
	ns[3, 3] = ns[1, 1] + ns[1, 2] + ns[2, 1] + ns[2, 2]
	
	#build return object and send it back
	return_obj = list()
	
	return_obj$results = raw_results
	return_obj$summary_square = summary_square
	return_obj$ns = ns
	return_obj$pct_data_used = round(ns[3, 3] / n, 3)
	return_obj$Ra = Ra
	return_obj$Rb = Rb
	return_obj$Rc = Rc
	return_obj$Rd = Rd
	sse = sum((raw_results$real_y - raw_results$est_true)^2, na.rm = TRUE)
	return_obj$oos_rmse = sqrt(sse / n)
	sst = sum((raw_results$real_y - mean(raw_results$real_y, na.rm = TRUE))^2, na.rm = TRUE)
	return_obj$out_of_sample_Rsq = 1 - sse / sst
	return_obj$pred_differences_avg = mean(abs(raw_results[, 1] - raw_results[, 2]), na.rm = TRUE)
	return_obj$pred_differences_sd = sd(abs(raw_results[, 1] - raw_results[, 2]), na.rm = TRUE)
	return_obj$avg_optimals = mean(c(Ra, Rd), na.rm = TRUE)
	return_obj$avg_non_optimals = mean(c(Rb, Rc), na.rm = TRUE)
	return_obj$avg_all = mean(c(Ra, Rb, Rc, Rd), na.rm = TRUE)
	avg_ys_tx_1 = mean(c(Ra, Rb), na.rm = TRUE)
	avg_ys_tx_2 = mean(c(Rc, Rd), na.rm = TRUE)
	if (avg_ys_tx_1 >= avg_ys_tx_2 && y_higher_is_better){ #sometimes continuous data aint continuous and you can have a "measure 0" event of equality here
		return_obj$avg_best = avg_ys_tx_1
	} else if (avg_ys_tx_1 >= avg_ys_tx_2 && !y_higher_is_better){
		return_obj$avg_best = avg_ys_tx_2
	} else if (avg_ys_tx_1 < avg_ys_tx_2 && y_higher_is_better){
		return_obj$avg_best = avg_ys_tx_2			
	} else if (avg_ys_tx_1 < avg_ys_tx_2 && !y_higher_is_better){
		return_obj$avg_best = avg_ys_tx_1
	}
	return_obj$q_adversarial = return_obj$avg_optimals - return_obj$avg_non_optimals
	return_obj$q_average = return_obj$avg_optimals - return_obj$avg_all
	return_obj$q_best = return_obj$avg_optimals - return_obj$avg_best
	return_obj	
}
