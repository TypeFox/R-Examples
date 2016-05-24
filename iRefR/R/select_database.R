####
# Select the subset of records in iRefIndex belonging to a certain primary Interaction Database:
####
select_database = function(database, MITAB_table, flag="this_database") {
	# 1. Dataset for these databases:
	dataset = NULL
	for (i in database) {
		dataset = rbind(dataset, MITAB_table[grep(i, MITAB_table$sourcedb, ignore.case=TRUE),])
	}

	# 2. Dataset for all except these databases:
	if (flag=="except_this_database") {
		if (dim(dataset)[1] > 0) {
			dataset = MITAB_table
			for (i in database) {
				dataset = dataset[-grep(i, dataset$sourcedb, ignore.case=TRUE),]
			}
		} else {
			dataset = MITAB_table
		}
	}

	output = dataset
}
