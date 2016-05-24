####
# Select the subset of records in iRefIndex belonging to a certain interaction type (binary, polymer,
# complex)
####
select_interaction_type = function(int_type, MITAB_table) {
	if (int_type == "binary") {
		output = MITAB_table[which(MITAB_table$edgetype == "X"),]
	}

	if (int_type == "polymer") {
		output = MITAB_table[which(MITAB_table$edgetype == "Y"),]
	}

	if (int_type == "complex") {
		output = MITAB_table[which(MITAB_table$edgetype == "C"),]
	}

	output
}
