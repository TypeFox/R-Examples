##################################################
# File:			rSCA.modeling.r
# Desp: 		R function for modeling using SCA
# Date: 		Jan 15, 2014, Regina, SK, Canada
# Author: 		Xiuquan Wang
# Email: 		xiuquan.wang@gmail.com
##################################################

##################################################
# Variables
# ================================================
#: new environment
rSCA.env = new.env()

#: significance level
rSCA.env$n_alpha = 0.95

#: data matrix
rSCA.env$o_sample_data_x = 0
rSCA.env$n_sample_size = 0
rSCA.env$n_sample_x_cols = 0
rSCA.env$o_sample_data_y = 0
rSCA.env$n_sample_y_cols = 0

#: matrix to store the sorted results
rSCA.env$o_sorted_matrix = 0
rSCA.env$o_sorted_temp_matrix = 0

#: output matrix
rSCA.env$o_output_tree = list()

#: stack to store the unprocessed node ids in the output tree
rSCA.env$o_nodeid_stack_cut = c()
rSCA.env$n_nodeid_statck_cut_cursor = 1
rSCA.env$o_nodeid_stack_merge = c()
rSCA.env$n_nodeid_statck_merge_cursor = 1

#: cutting and merging flags for loop, 1:do loop, 0:no need
rSCA.env$n_flag_cut = 0
rSCA.env$n_flag_merge = 0

#: some statistical infos for the results
rSCA.env$n_cut_times = 0
rSCA.env$n_merge_times = 0
rSCA.env$n_leafnodes_count = 0

#: hash table, hash function = (a + b) % n_sample_size
#: structure: [value, a, b, pointer]
rSCA.env$o_hashtable_matrix = matrix(0, rSCA.env$n_sample_size, 4)
#: hash table list: (value, a, b, pointer)
rSCA.env$o_hashtable_list = list()
rSCA.env$n_hashtable_list_index = 0

#: output files & paths
rSCA.env$s_tree_file = ""
rSCA.env$s_tree_filepath = ""
rSCA.env$s_map_file = ""
rSCA.env$s_map_filepath = ""
rSCA.env$s_logfile = ""
rSCA.env$s_logfilepath = ""

#: map representation: interval or mean
rSCA.env$n_mapvalue = "mean"

#: using optimization seeking method (Golden Section Search) or not
rSCA.env$b_GSS = FALSE

#: debug display
rSCA.env$b_debug = FALSE

##################################################
#
# NAME:						rSCA.modeling
#
# INPUTs:
#
#		@alpha:				significance level, usually in 0.001 - 0.05, default is 0.05
#
#		@xfile: 			a string to specify the full file name of the x file, only supports *.txt or *.csv
#
#		@yfile: 			for y file, similar to x file
#
#		@x.row.names:		TRUE/FALSE, default is FALSE
#
#		@x.col.names:		TRUE/FALSE,	default is FALSE
#
#		@y.row.names:		TRUE/FALSE, default is FALSE
#
#		@y.col.names:		TRUE/FALSE,	default is FALSE
#
#		@x.missing.flag:	a string to specify the missing flag, default is "NA"
#
#		@y.missing.flag:	a string to specify the missing flag, default is "NA"
#
#		@x.type:			".txt" or ".csv", default is ".txt"
#
#		@y.type:			".txt" or ".csv", default is ".txt"
#
#		@mapvalue:			{mean, max, min, median, interval, radius, variation, random}
#							default is mean, 
#							interval: [min, max], 
#							radius: [mean, radius], radius = (max - min) / 2,
#							variation: [mean, sd], 
#							random: generate randomly from [min, max].
#
#		@GSS:				TRUE/FALSE, default is FALSE (no optimization seeking), the Golden Section Search method is used for seeking the best cutting point if GSS is TRUE
#
#		@debug:				TRUE/FALSE, default is FALSE, log file will be created if debug is TRUE
#
# OUTPUTs:
#
#		@statistics infos:	to screen
#
#		@detailed infos:	to log file, if debug is TRUE
#
# RETURNS:
#
#		@model object:		treefile + mapfile + type (mean/interval)
#
###################################################
rSCA.modeling <- function(alpha = 0.05, xfile, yfile, x.row.names = FALSE, x.col.names = FALSE, y.row.names = FALSE, y.col.names = FALSE, x.missing.flag = "NA", y.missing.flag = "NA", x.type = ".txt", y.type = ".txt", mapvalue = "mean", GSS = FALSE, debug = FALSE)
{
	rSCA.env$n_alpha = alpha
	rSCA.env$n_mapvalue = mapvalue
	rSCA.env$b_GSS = GSS
	rSCA.env$b_debug = debug

	#: model list
	o_model = list(treefile = "", mapfile = "", logfile = "", type = "mean", totalNodes = 0, leafNodes = 0, cuttingActions = 0, mergingActions = 0)

	#: data matrix
	o_xdata = 0
	o_ydata = 0
	
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
	
	#: read y data file
	if (y.type == ".txt")
	{
		if (y.row.names == TRUE && y.col.names == TRUE)
			o_ydata = read.table(yfile, header = TRUE, row.names = 1, na.strings = y.missing.flag)
		else if (y.row.names == TRUE && y.col.names == FALSE)
			o_ydata = read.table(yfile, header = FALSE, row.names = 1, na.strings = y.missing.flag)
		else if (y.row.names == FALSE && y.col.names == TRUE)
			o_ydata = read.table(yfile, header = TRUE, na.strings = y.missing.flag)
		else if (y.row.names == FALSE && y.col.names == FALSE)
			o_ydata = read.table(yfile, header = FALSE, na.strings = y.missing.flag)
	}
	if (y.type == ".csv")
	{
		if (y.row.names == TRUE && y.col.names == TRUE)
			o_ydata = read.csv(yfile, header = TRUE, row.names = 1, na.strings = y.missing.flag)
		else if (y.row.names == TRUE && y.col.names == FALSE)
			o_ydata = read.csv(yfile, header = FALSE, row.names = 1, na.strings = y.missing.flag)
		else if (y.row.names == FALSE && y.col.names == TRUE)
			o_ydata = read.csv(yfile, header = TRUE, na.strings = y.missing.flag)
		else if (y.row.names == FALSE && y.col.names == FALSE)
			o_ydata = read.csv(yfile, header = FALSE, na.strings = y.missing.flag)
	}
	
	#: remove missing rows
	o_xdata = na.omit(o_xdata)
	o_ydata = na.omit(o_ydata)
	
	#: basic infos
	n_xdatarow = nrow(o_xdata)
	n_ydatarow = nrow(o_ydata)
	
	#: checking if x and y have the same number of rows.
	if (n_xdatarow != n_ydatarow)
	{
		cat("Error: please make sure x and y have the same number of rows.\n")
	}
	else
	{
		#: start training
		rSCA.env$o_sample_data_x = o_xdata
		rSCA.env$n_sample_size = nrow(rSCA.env$o_sample_data_x)
		rSCA.env$n_sample_x_cols = ncol(rSCA.env$o_sample_data_x)

		rSCA.env$o_sample_data_y = o_ydata
		rSCA.env$n_sample_y_cols = ncol(rSCA.env$o_sample_data_y)
		
		s_tmp_filename = get_random_filename()
		s_tmp_treefile = paste("tree_", s_tmp_filename, ".txt", sep = "")
		s_tmp_mapfile = paste("map_", s_tmp_filename, ".txt", sep = "")
		s_tmp_logfile = paste("log_", s_tmp_filename, ".txt", sep = "")

		rSCA.env$s_tree_file = paste("tree_", s_tmp_filename, ".txt", sep = "")
		rSCA.env$s_tree_filepath = rSCA.env$s_tree_file
		rSCA.env$s_map_file = paste("map_", s_tmp_filename, ".txt", sep = "")
		rSCA.env$s_map_filepath = rSCA.env$s_map_file
		rSCA.env$s_logfile = paste("log_", s_tmp_filename, ".txt", sep = "")
		rSCA.env$s_logfilepath = rSCA.env$s_logfile

		#: start SCA modeling
		do_cluster()
		
		#: return the model
		o_model$treefile = rSCA.env$s_tree_file
		o_model$mapfile = rSCA.env$s_map_file
		if (rSCA.env$b_debug)
			o_model$logfile = rSCA.env$s_logfile
		else
			o_model$logfile = "NA"
		o_model$type = rSCA.env$n_mapvalue
		o_model$totalNodes = length(rSCA.env$o_output_tree)
		o_model$leafNodes = rSCA.env$n_leafnodes_count
		o_model$cuttingActions = rSCA.env$n_cut_times
		o_model$mergingActions = rSCA.env$n_merge_times
		
		return(o_model)
	}
}


