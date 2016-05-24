##################################################
# File:			rSCA.missing.r
# Desp: 		R function for missing check
# Date: 		Jan 14, 2014, Regina, SK, Canada
# Author: 		Xiuquan Wang
# Email: 		xiuquan.wang@gmail.com
##################################################


##################################################
#
# NAME:					rSCA.missing
#
# INPUTs:
#
#		@file: 			a string to specify the full file name of the data file, only supports *.txt or *.csv
#
#		@row.names:		TRUE/FALSE, default is FALSE
#
#		@col.names:		TRUE/FALSE,	default is FALSE
#
#		@missing.flag:	a string to specify the missing flag, default is "NA"
#
#		@type:			".txt" or ".csv", default is ".txt"
#
# OUTPUTs:
#
#		@statistics:	 				missing		total		percent
#						by element:		10			100			10 %
#						by row:			5			10			50 %
#						by column:		5			10			50 %
#
#		@location:		     [,1] [,2] [,3] [,4]
#						[1,]  -    -    -    -  
#						[2,]  -   <?>   -    -  
#						[3,]  -    -    -    -  
#
#
# RETURNS:
#
#		@bool:			TRUE or FALSE, indicates if it passes the missing check
#
###################################################
rSCA.missing <- function(file, row.names = FALSE, col.names = FALSE, missing.flag = "NA", type = ".txt")
{
	#: pass flag
	b_pass = FALSE
	
	#: data matrix
	o_data = 0
	
	#: read data file
	if (type == ".txt")
	{
		if (row.names == TRUE && col.names == TRUE)
			o_data = read.table(file, header = TRUE, row.names = 1, na.strings = missing.flag)
		else if (row.names == TRUE && col.names == FALSE)
			o_data = read.table(file, header = FALSE, row.names = 1, na.strings = missing.flag)
		else if (row.names == FALSE && col.names == TRUE)
			o_data = read.table(file, header = TRUE, na.strings = missing.flag)
		else if (row.names == FALSE && col.names == FALSE)
			o_data = read.table(file, header = FALSE, na.strings = missing.flag)
	}
	if (type == ".csv")
	{
		if (row.names == TRUE && col.names == TRUE)
			o_data = read.csv(file, header = TRUE, row.names = 1, na.strings = missing.flag)
		else if (row.names == TRUE && col.names == FALSE)
			o_data = read.csv(file, header = FALSE, row.names = 1, na.strings = missing.flag)
		else if (row.names == FALSE && col.names == TRUE)
			o_data = read.csv(file, header = TRUE, na.strings = missing.flag)
		else if (row.names == FALSE && col.names == FALSE)
			o_data = read.csv(file, header = FALSE, na.strings = missing.flag)
	}
	
	#: basic infos
	n_datarow = nrow(o_data)
	n_datacol = ncol(o_data)
	n_totaldata = n_datarow * n_datacol
	
	#: convert to True or False (NA --> True)
	o_missing_matrix = is.na(o_data)
	n_totalmissing = sum(o_missing_matrix)
	n_totalmissing_percent = round((n_totalmissing / n_totaldata) * 100, digits = 2)
	
	#: calculate missing rows
	n_row_missing = 0
	n_row_missing_percent = 0
	for (irow in 1:n_datarow)
	{
		if (sum(o_missing_matrix[irow, ]) > 0)
			n_row_missing = n_row_missing + 1
	}
	n_row_missing_percent = round((n_row_missing / n_datarow) * 100, digits = 2)
	
	#: calculate missing cols
	n_col_missing = 0
	n_col_missing_percent = 0
	for (icol in 1:n_datacol)
	{
		if (sum(o_missing_matrix[ , icol]) > 0)
			n_col_missing = n_col_missing + 1
	}
	n_col_missing_percent = round((n_col_missing / n_datacol) * 100, digits = 2)
	
	#: construct missing statistics
	s_missing_statistics = matrix(0, 3, 3)
	rownames(s_missing_statistics) = c("by element: ", "by row: ", "by column: ")
	colnames(s_missing_statistics) = c("total ", "missing ", "percent")
	s_missing_statistics[1, ] = c(n_totaldata, n_totalmissing, paste(n_totalmissing_percent, "%", sep = " "))
	s_missing_statistics[2, ] = c(n_datarow, n_row_missing, paste(n_row_missing_percent, "%", sep = " "))
	s_missing_statistics[3, ] = c(n_datacol, n_col_missing, paste(n_col_missing_percent, "%", sep = " "))
	
	#: construct missing table
	s_missing_table = matrix(" - ", n_datarow, n_datacol)
	for (i in 1:n_datarow)
	{
		for (j in 1:n_datacol)
		{
			if (o_missing_matrix[i, j])
			{
				s_missing_table[i, j] = "<?>"
			}
		}
	}
	
	#: print outputs
	cat("\nStatistics:\n\n")
	print(s_missing_statistics, quote = FALSE, zero.print = ".")
	cat("\n\nMissing table:\n\n")
	print(s_missing_table, quote = FALSE)
	cat("\n")
	
	#: check if it got passed
	if (n_totalmissing <= 0)	b_pass = TRUE
	
	#: return
	return(b_pass)
}
