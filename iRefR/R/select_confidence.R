####
# Select the subset of records of iRefIndex corresponding to a given set of confidence scores:
####
select_confidence = function(confidence_type="lpr", confidence_value, MITAB_table) {
	# 1. Generate a vector of the strings to grep:
	vector_to_grep = NULL
	for (i in confidence_value) {
		vector_to_grep = c(vector_to_grep, paste(confidence_type, "\\:", i, "\\|", sep=""))
	}

	# 2. Get MITAB rows with specified confidence:
	MITAB_list = list()
	for (i in vector_to_grep) {
		MITAB_list[[i]] = MITAB_table[grep(i, MITAB_table$confidence),]
	}

	MITAB_output = NULL
	for (i in vector_to_grep) {
		MITAB_output = rbind(MITAB_output, MITAB_list[[i]])
	}

	MITAB_output
}
