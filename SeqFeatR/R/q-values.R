#R Script written by Bettina Budeus, rokkaku-fight@gmx.de
#from 2012-07-23 to 2012--
#calculates and adds q_values to a list of given p_values in csv format with p_value as column header for this column.

#search for "change" for points where changes are possible
#please check before executing

#libs
library(qvalue)

q_value_set_input_file <- structure(function(
	### set the input file for calculating the q-values
	path_to_file
	### file with a column with p-values. Has to be in csv and contain a column called p_value.  
	### For reference please look in example file.	
	){
	data <- read.csv2(path_to_file,na.strings = "", stringsAsFactors = FALSE)
	.GlobalEnv[["input_matrix"]] <- data
},ex=function(){
	q_value_set_input_file("SeqFeatR/extdata/co_mutation_results.csv")
})

qvalues <- structure(function(# Calculate q-values
	### calculates and adds q_values to a list of given p_values in csv format with p_value as column header for this column.
	##details<< Takes a csv file with a column which is called p_values and uses the qvalues package to calculate from this column the corresponding q-values.
	## Uses qvalue package and the calculation within this package to estimate the q-values.
	path_to_file_csv = NULL,
	### file with a column with p-values. Has to be in csv and contain a column called p_value.  
	### For reference please look in example file.	
	save_name_csv
	### output file name	
	){
	result <- calculate_qvalues_inner(path_to_file_csv, save_name_csv)
	return (result)
},ex=function(){
	mut <- system.file("extdata", "co_mutation_results.csv", package="SeqFeatR")
	qvalues(mut, "csv_with_q_values.csv")
})

calculate_qvalues_inner <- function(path_to_file = NULL, save_name){
	if (is.null(path_to_file)==FALSE){
		q_value_set_input_file(path_to_file)
	}

	input_matrix <- .GlobalEnv[["input_matrix"]]
	
	### contains the p_values of the loaded file
	p_values <- as.numeric(input_matrix$p_value)
	print (min(p_values))
	print (max(p_values))
	### contains the absolute p_values of the loaded file
	p_v <- c()

	for (entry in 1:length(p_values)){
		### the absolut p_value
		value <- abs(as.numeric(p_values[entry]))
		p_v <- c(p_v,value)
	}
	### the calculatet list of p_values
	q_value <- qvalue(p=p_v, lambda=0, fdr.level=0.05)
	### the result of the analysis with added column for q-values
	result_matrix <- cbind(input_matrix, q_value$qvalues)

	colnames(result_matrix)[ncol(result_matrix)] <- "q_value"
	
	result_matrix <- result_matrix[order(result_matrix$q_value),]

	write.csv2(result_matrix, paste(save_name, sep=""))	

	return (result_matrix)
	### returns the input file as matrix with added column with q-values
}

