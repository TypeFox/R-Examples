##################################################
# File:			rSCA.inference.r
# Desp: 		R function for inference with SCA
# Date: 		Jan 16, 2014, Regina, SK, Canada
# Author: 		Xiuquan Wang
# Email: 		xiuquan.wang@gmail.com
##################################################

#: load rSCA.prediction.r
#source("rSCA.prediction.r")

##################################################
# Variables
# ================================================
#: new environment
rSCA.env = new.env()

rSCA.env$o_result_tree_p = 0  #o_result_tree = 0
rSCA.env$n_result_tree_rows_p = 0 #n_result_tree_rows_p = 0
rSCA.env$o_mean_y_p = 0
rSCA.env$n_y_cols_p = 0
rSCA.env$o_predictors_p = 0
rSCA.env$n_predictors_rows_p = 0
rSCA.env$o_predictants_p = 0

rSCA.env$s_result_file_p = ""
rSCA.env$s_result_filepath_p = ""

rSCA.env$n_model_type_p = ""

##################################################
#
# NAME:						rSCA.inference
#
# INPUTs:
#
#		@xfile: 			a string to specify the full file name of the x file, only supports *.txt or *.csv
#
#		@x.row.names:		TRUE/FALSE, default is FALSE
#
#		@x.col.names:		TRUE/FALSE, default is FALSE
#
#		@x.missing.flag:	a string to specify the missing flag, default is "NA"
#
#		@x.type:			".txt" or ".csv", default is ".txt"
#
#		@model:				a list object pointing to the SCA model
#
# OUTPUTs:
#
#		@result file:		a text file contains the predicted values
#
# RETURNS:
#
#		no return.
#
###################################################
rSCA.inference <- function(xfile, x.row.names = FALSE, x.col.names = FALSE, x.missing.flag = "NA", x.type = ".txt", model)
{
	#: data matrix
	o_xdata = 0
	
	#: read x data file
	if (x.type == ".txt")
	{
		if (x.row.names == TRUE && x.col.names == TRUE)
			o_xdata = read.table(xfile, header = TRUE, row.names = 1, na.strings = x.missing.flag)
		else if (x.row.names == TRUE && x.col.names == FALSE)
			o_xdata = read.table(xfile, header = FALSE, row.names = 1, na.strings = x.missing.flag)
		else if (x.row.names == FALSE && x.col.names == TRUE)
			o_xdata = read.table(xfile, header = TRUE, na.strings = x.missing.flag)
		else if (x.row.names == FALSE && x.col.names == FALSE)
			o_xdata = read.table(xfile, header = FALSE, na.strings = x.missing.flag)
	}
	if (x.type == ".csv")
	{
		if (x.row.names == TRUE && x.col.names == TRUE)
			o_xdata = read.csv(xfile, header = TRUE, row.names = 1, na.strings = x.missing.flag)
		else if (x.row.names == TRUE && x.col.names == FALSE)
			o_xdata = read.csv(xfile, header = FALSE, row.names = 1, na.strings = x.missing.flag)
		else if (x.row.names == FALSE && x.col.names == TRUE)
			o_xdata = read.csv(xfile, header = TRUE, na.strings = x.missing.flag)
		else if (x.row.names == FALSE && x.col.names == FALSE)
			o_xdata = read.csv(xfile, header = FALSE, na.strings = x.missing.flag)
	}

	#: remove missing rows
	o_xdata = na.omit(o_xdata)
	rSCA.env$o_predictors_p = o_xdata
	rSCA.env$n_predictors_rows_p = nrow(o_xdata)
	
	#: read tree and map file
	rSCA.env$o_result_tree_p = read.table(model$treefile, header = TRUE)
	rSCA.env$n_result_tree_rows_p = nrow(rSCA.env$o_result_tree_p)
	rSCA.env$o_mean_y_p = read.table(model$mapfile, header = TRUE)
	rSCA.env$n_y_cols_p = ncol(rSCA.env$o_mean_y_p)
	
	#: define the result file
	s_tmp_file = paste("rsl_", get_random_filename(), ".txt", sep = "")
	rSCA.env$s_result_file_p = s_tmp_file
	rSCA.env$s_result_filepath_p = s_tmp_file
	
	#: set model type
	rSCA.env$n_model_type_p = model$type

	#: start prediction
	do_prediction()

}
