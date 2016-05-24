print.PTE_bootstrap_results = function(x, ...){
	cat("    I_adversarial observed est = ", 
			round(x$observed_q_adversarial, 3), 
			",  p val = ",
			round(x$p_val_adversarial, 3), 
			", \n      ",			
			round(100 * (1 - x$alpha), 1),
			"% CI's: pctile = [",
			round(x$ci_q_adversarial[1], 3), 
			", ",
			round(x$ci_q_adversarial[2], 3), 
			"], ",
			"BCa = [",
			round(x$bca_ci_q_adversarial[1], 3), 
			", ",
			round(x$bca_ci_q_adversarial[2], 3), 
			"],", 
			sep = "")
	cat("\n    I_random observed_est = ", 
			round(x$observed_q_average, 3),
			",  p val = ", 
			round(x$p_val_average, 3),
			", \n      ",
			round(100 * (1 - x$alpha), 1),
			"% CI's: pctile = [",
			round(x$ci_q_average[1], 3), 
			", ",
			round(x$ci_q_average[2], 3), 
			"], ",
			"BCa = [",
			round(x$bca_ci_q_average[1], 3), 
			", ",
			round(x$bca_ci_q_average[2], 3), 
			"],",
			sep = "")
	cat("\n    I_best observed_est = ", 
			round(x$observed_q_best, 3), 
			",  p val = ",
			round(x$p_val_best, 3),
			", \n      ",
			round(100 * (1 - x$alpha), 1),
			"% CI's: pctile = [",
			round(x$ci_q_best[1], 3), 
			", ",
			round(x$ci_q_best[2], 3),  
			"], ",
			"BCa = [",
			round(x$bca_ci_q_best[1], 3), 
			", ",
			round(x$bca_ci_q_best[2], 3),  
			"]", 
			sep = "")
	cat("\n")  
}

#alias the print function
summary.PTE_bootstrap_results = function(object, ...){
	print(object)
}