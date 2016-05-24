##################################################
# File:			rSCA.correlation.r
# Desp: 		R function for correlation check
# Date: 		Jan 15, 2014, Regina, SK, Canada
# Author: 		Xiuquan Wang
# Email: 		xiuquan.wang@gmail.com
##################################################


##################################################
#
# NAME:						rSCA.correlation
#
# INPUTs:
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
# OUTPUTs:
#
#		@correlation table:			  	x.V1 	x.V2	x.V3 	x.V4
#								y.V1  	0.8 	0.9		0.78	0.23
#								y.V2  	0.9 	0.5		0.48	0.7
#								y.V3  	0.3 	0.2		0.3		0.63
#
# RETURNS:
#
#		no return.
#
###################################################
rSCA.correlation <- function(xfile, yfile, x.row.names = FALSE, x.col.names = FALSE, y.row.names = FALSE, y.col.names = FALSE, x.missing.flag = "NA", y.missing.flag = "NA", x.type = ".txt", y.type = ".txt")
{
	
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
		#: contruct correlation table
		o_correlation_table = round(cor(o_ydata, o_xdata), digits = 4)
		rownames(o_correlation_table) = colnames(o_ydata)
		colnames(o_correlation_table) = colnames(o_xdata)
		
		cat("\nCorrelation table:\n\n")
		print(o_correlation_table)
		cat("\n")
	}
}
