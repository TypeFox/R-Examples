##########################
# Merge tables in complexList format:
##########################
merge_complexes_lists = function(list_of_complexLists) {
	# 1. Making unique tables for each complexList:
	n = 0; new_table = list(); new_list_of_complexLists = list()
	for (i in list_of_complexLists) {
		n = n + 1
		colnames(i) = NULL
		list_unique_complexes = unique(i[,2])
		for (j in 1:length(list_unique_complexes)) {
			tmp = i[which(i[,2]==list_unique_complexes[j]),]
			if (is.matrix(tmp)=="FALSE") {
				new_table[[j]] = tmp[1:2]
			} else {
				new_table[[j]] = tmp[1,]
			}
		}
		new_list_of_complexLists[[n]] = do.call(rbind, new_table)
	}

	# 2. Function to remove from less prioritized list:
	remotionator = function(list1, list2) {
		matches = which((list2[,2] %in% list1[,2])==TRUE)
		if (length(matches)>0) {
			new_list2 = list2[-matches,]
		} else {
			new_list2 = list2
		}
		new_list2
	}

	# 3. Actual comparison and prioritized remotion of repeated complexes:
	n = length(new_list_of_complexLists)
	for (i in 1:(n-1)) {
		for (j in (i+1):n) {
			new_list_of_complexLists[[j]] = remotionator(new_list_of_complexLists[[i]], new_list_of_complexLists[[j]])
		}
	}

	merged_table = do.call(rbind, new_list_of_complexLists)
	colnames(merged_table) = c("complex ID", "subunits")

	result = merged_table

}
