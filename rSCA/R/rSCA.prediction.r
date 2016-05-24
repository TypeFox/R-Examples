#################################################################
# Filename: 	rSCA.prediction.r
# Date: 		2010/11/06, Regina, SK, Canada
# Author: 		Xiuquan Wang
# Email:		xiuquan.wang@gmail.com
# Desp: 		Doing prediction using the SCA model
# ===============================================================
# History: 		2010/11/06	created by Xiuquan Wang
##################################################################

# ---------------------------------------------------------------
# Function definitions
# ---------------------------------------------------------------

#: predict one input
#: input: o_precitor = c(x1, x2, x3, ... )
#: output: o_predictant = c(y1, y2, y3, ... )
f_predict_one = function(o_precitor)
{
	o_predictant = c()

	#: initiate as c(0, 0, 0, ... )
	for (ip in 1:rSCA.env$n_y_cols_p)
	{
		o_predictant[ip] = 0
	}

	n_tree_index = 1
	while(n_tree_index <= rSCA.env$n_result_tree_rows_p)
	{
		if (rSCA.env$o_result_tree_p[n_tree_index, 4] == -1 && rSCA.env$o_result_tree_p[n_tree_index, 5] == -1)
		{
			#: if meets leaf, then get the mean y
			o_predictant = c(as.numeric(rSCA.env$o_mean_y_p[rSCA.env$o_result_tree_p[n_tree_index, 1], ]))
			break
		}

		#: if not, compare with the Xij
		#: if <= it, go to left node; otherwise go to right node

		#: get j (column index)
		n_j = rSCA.env$o_result_tree_p[n_tree_index, 2]
		if (n_j <= 0)
		{
			#: if merged, just go to the left node
			n_tree_index = rSCA.env$o_result_tree_p[n_tree_index, 4]
			next
		}
		if (o_precitor[n_j] <= rSCA.env$o_result_tree_p[n_tree_index, 3])
		{
			n_tree_index = rSCA.env$o_result_tree_p[n_tree_index, 4]
		}
		else
		{
			n_tree_index = rSCA.env$o_result_tree_p[n_tree_index, 5]
		}
	}
	
	#: return
	return(o_predictant)
}

#: do precition to all input
f_predict = function()
{
	for (iPredictor in 1:rSCA.env$n_predictors_rows_p)
	{
		o_vector_predictor = c(as.numeric(rSCA.env$o_predictors_p[iPredictor, ]))
		o_vector_predictant = f_predict_one(o_vector_predictor)
		
		rSCA.env$o_predictants_p[iPredictor, ] <- o_vector_predictant
	}
}

# ---------------------------------------------------------------
# Initialization function
# ---------------------------------------------------------------
f_init_p = function()
{
	cat("Initializing...\t\t")

	#: predictant to be outputed
	rSCA.env$o_predictants_p <- matrix(0, rSCA.env$n_predictors_rows_p, rSCA.env$n_y_cols_p)

	cat("SUCCESS!\r\n")
}

# ---------------------------------------------------------------
# Main funtion
# ---------------------------------------------------------------
f_main_p = function()
{
	cat("Processing...\t\t")
	f_predict()
	cat("SUCCESS!\r\n")

	#: set column names + format output
	if (rSCA.env$n_model_type_p == "mean" || rSCA.env$n_model_type_p == "max" || rSCA.env$n_model_type_p == "min" || rSCA.env$n_model_type_p == "median")
	{
		colnames(rSCA.env$o_predictants_p) = colnames(rSCA.env$o_mean_y_p)
		write.table(rSCA.env$o_predictants_p, file = rSCA.env$s_result_filepath_p, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
	}
	if (rSCA.env$n_model_type_p == "interval")
	{
		#: print col names
		s_colnames = colnames(rSCA.env$o_mean_y_p)
		n_midcol = rSCA.env$n_y_cols_p / 2
		for (icn in 1:n_midcol)
		{
			cat(s_colnames[icn], "\t\t", file = rSCA.env$s_result_filepath_p, sep = "", append = TRUE)
		}
		cat("\r\n", file = rSCA.env$s_result_filepath_p, sep = "", append = TRUE)
		
		#: print result matrix
		for (ime in 1:nrow(rSCA.env$o_predictants_p))
		{
			o_vec_min = rSCA.env$o_predictants_p[ime, 1:n_midcol]
			o_vec_max = rSCA.env$o_predictants_p[ime, (n_midcol + 1):rSCA.env$n_y_cols_p]
			for (iv in 1:n_midcol)
			{
				cat("[", o_vec_min[iv], ", ", o_vec_max[iv], "]\t", file = rSCA.env$s_result_filepath_p, sep = "", append = TRUE)
			}
			cat("\r\n", file = rSCA.env$s_result_filepath_p, sep = "", append = TRUE)
		}
	}
	if (rSCA.env$n_model_type_p == "radius")
	{
		#: print col names
		s_colnames = colnames(rSCA.env$o_mean_y_p)
		n_midcol = rSCA.env$n_y_cols_p / 2
		for (icn in 1:n_midcol)
		{
			cat(s_colnames[icn], "\t\t", file = rSCA.env$s_result_filepath_p, sep = "", append = TRUE)
		}
		cat("\r\n", file = rSCA.env$s_result_filepath_p, sep = "", append = TRUE)
		
		#: print result matrix
		for (ime in 1:nrow(rSCA.env$o_predictants_p))
		{
			o_vec_mean = rSCA.env$o_predictants_p[ime, 1:n_midcol]
			o_vec_radius = rSCA.env$o_predictants_p[ime, (n_midcol + 1):rSCA.env$n_y_cols_p]
			for (iv in 1:n_midcol)
			{
				cat("[", o_vec_mean[iv], " +/- ", o_vec_radius[iv], "]\t", file = rSCA.env$s_result_filepath_p, sep = "", append = TRUE)
			}
			cat("\r\n", file = rSCA.env$s_result_filepath_p, sep = "", append = TRUE)
		}
	}
	if (rSCA.env$n_model_type_p == "variation")
	{
		#: print col names
		s_colnames = colnames(rSCA.env$o_mean_y_p)
		n_midcol = rSCA.env$n_y_cols_p / 2
		for (icn in 1:n_midcol)
		{
			cat(s_colnames[icn], "\t\t", file = rSCA.env$s_result_filepath_p, sep = "", append = TRUE)
		}
		cat("\r\n", file = rSCA.env$s_result_filepath_p, sep = "", append = TRUE)
		
		#: print result matrix
		for (ime in 1:nrow(rSCA.env$o_predictants_p))
		{
			o_vec_mean = rSCA.env$o_predictants_p[ime, 1:n_midcol]
			o_vec_sd = rSCA.env$o_predictants_p[ime, (n_midcol + 1):rSCA.env$n_y_cols_p]
			for (iv in 1:n_midcol)
			{
				cat("[", o_vec_mean[iv], " +/- ", o_vec_sd[iv], "]\t", file = rSCA.env$s_result_filepath_p, sep = "", append = TRUE)
			}
			cat("\r\n", file = rSCA.env$s_result_filepath_p, sep = "", append = TRUE)
		}
	}
	if (rSCA.env$n_model_type_p == "random")
	{
		#: print col names
		s_colnames = colnames(rSCA.env$o_mean_y_p)
		n_midcol = rSCA.env$n_y_cols_p / 2
		for (icn in 1:n_midcol)
		{
			cat(s_colnames[icn], "\t\t", file = rSCA.env$s_result_filepath_p, sep = "", append = TRUE)
		}
		cat("\r\n", file = rSCA.env$s_result_filepath_p, sep = "", append = TRUE)
		
		#: print result matrix
		for (ime in 1:nrow(rSCA.env$o_predictants_p))
		{
			o_vec_min = rSCA.env$o_predictants_p[ime, 1:n_midcol]
			o_vec_max = rSCA.env$o_predictants_p[ime, (n_midcol + 1):rSCA.env$n_y_cols_p]
			for (iv in 1:n_midcol)
			{
				n_random_value = runif(1, o_vec_min[iv], o_vec_max[iv])
				cat(n_random_value, "\t", file = rSCA.env$s_result_filepath_p, sep = "", append = TRUE)
			}
			cat("\r\n", file = rSCA.env$s_result_filepath_p, sep = "", append = TRUE)
		}
	}
	
	#: print to screen
	cat("Result File:\t\t", rSCA.env$s_result_filepath_p, "\r\n", sep = "")
}

# ---------------------------------------------------------------
# Interface function
# ---------------------------------------------------------------
do_prediction = function()
{
	# : store the start time
	time_stat <- proc.time()

	#: initialize
	f_init_p()

	#: do main function
	f_main_p()

	# : calculate the total time used
	time_end <- (proc.time() - time_stat)[[3]]
	Hours <- time_end %/% (60*60)
	Minutes <- (time_end %% 3600) %/% 60
	Seconds <- time_end %% 60
	time_used <- paste(Hours, " h ", Minutes, " m ", Seconds, " s.", sep="")
	cat("Time Used:\t\t", time_used, "\r\n", sep = "")
}


