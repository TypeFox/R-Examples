####
# Convert a table from edgeList format to MITAB format:
####
convert_edgeList_to_MITAB = function(edgeList, MITAB_table, canonical_rep="yes", int_type="all") {

	# 1. Separate complexes from binaries and polymers:
	MITAB_binary = select_interaction_type("binary", MITAB_table)
	MITAB_complex = select_interaction_type("complex", MITAB_table)
	MITAB_polymer = select_interaction_type("polymer", MITAB_table)
	possible_binaries = NULL
	possible_complexes = NULL
	possible_polymers = NULL

	# 2. Reconstructing binaries:
	if (int_type=="binary" || int_type=="all") {
		if (canonical_rep == "no") {
			possible_binary_1 = MITAB_binary[which(MITAB_binary$irogida %in% edgeList[,1] & MITAB_binary$irogidb %in% edgeList[,2]),]
			possible_binary_2 = MITAB_binary[which(MITAB_binary$irogidb %in% edgeList[,1] & MITAB_binary$irogida %in% edgeList[,2]),]
			possible_binaries = unique(rbind(possible_binary_1, possible_binary_2))
		} else {
			possible_binary_1 = MITAB_binary[which(MITAB_binary$icrogida %in% edgeList[,1] & MITAB_binary$icrogidb %in% edgeList[,2]),]
			possible_binary_2 = MITAB_binary[which(MITAB_binary$icrogidb %in% edgeList[,1] & MITAB_binary$icrogida %in% edgeList[,2]),]
			possible_binaries = unique(rbind(possible_binary_1, possible_binary_2))
		}
	}

	# 3. Reconstructing complexes:
	if (int_type=="complex" || int_type=="all") {
		rigid_list = NULL
		for (i in 1:dim(edgeList)[1]) {
			if (canonical_rep == "no") {
				rigid_node_a = MITAB_complex[which(edgeList[i,1]==MITAB_complex$irogidb), "irigid"]
				rigid_node_b = MITAB_complex[which(edgeList[i,2]==MITAB_complex$irogidb), "irigid"]
			} else {
				rigid_node_a = MITAB_complex[which(edgeList[i,1]==MITAB_complex$icrogidb), "icrigid"]
				rigid_node_b = MITAB_complex[which(edgeList[i,2]==MITAB_complex$icrogidb), "icrigid"]
			}
			intersection = intersect(rigid_node_a, rigid_node_b)
			if (length(intersection) > 0) {
				rigid_list = c(rigid_list, intersection)
			}
		}

		if (canonical_rep == "no") {
			possible_complexes = MITAB_complex[which(MITAB_complex$irigid %in% rigid_list),]
		} else {
			possible_complexes = MITAB_complex[which(MITAB_complex$icrigid %in% rigid_list),]
		}
	}

	# 4. Reconstructing possible polymers:
	if (int_type=="polymer" || int_type=="all") {
		if (canonical_rep == "no") {
			possible_polymers = MITAB_polymer[which(MITAB_polymer$irogida %in% edgeList[,1] & MITAB_polymer$irogidb %in% edgeList[,2]),]
		} else {
			possible_polymers = MITAB_polymer[which(MITAB_polymer$icrogida %in% edgeList[,1] & MITAB_polymer$icrogidb %in% edgeList[,2]),]
		}
	}

	output = rbind(possible_binaries, possible_complexes, possible_polymers)
}
