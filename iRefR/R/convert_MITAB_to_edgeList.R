####
# Convert a table from MITAB format to edgeList format:
####
convert_MITAB_to_edgeList = function(MITAB_table, edge_weight_col="default", complex_rep="spoke", canonical_rep="yes", directionality="undirected", node_names="rogs", multi_edge="no") {
	# 1. Separate complexes from binaries and polymers:
	if (edge_weight_col=="default") {
		edge_weight_col = 55
		MITAB_table[,edge_weight_col] = rep(1, dim(MITAB_table)[1])
	}
	MITAB_binary = select_interaction_type("binary", MITAB_table)
	MITAB_complex = select_interaction_type("complex", MITAB_table)
	MITAB_polymer = select_interaction_type("polymer", MITAB_table)

	# 2. Getting uid, irogid or icrogid columns:
	if (node_names == "uids") {
		edges_binary = unique(cbind(as.character(MITAB_binary$X.uidA), as.character(MITAB_binary$uidB), as.character(MITAB_binary$experimental_role_A), as.character(MITAB_binary$experimental_role_B), MITAB_binary$irigid, as.character(MITAB_binary$edgetype), MITAB_binary[,edge_weight_col]))
		edges_complex = unique(cbind(as.character(MITAB_complex$X.uidA), as.character(MITAB_complex$uidB), as.character(MITAB_complex$experimental_role_A), as.character(MITAB_complex$experimental_role_B), MITAB_complex$irigid, as.character(MITAB_complex$edgetype), MITAB_complex[,edge_weight_col]))
		edges_polymer = unique(cbind(as.character(MITAB_polymer$X.uidA), as.character(MITAB_polymer$uidB), MITAB_polymer$irigid, as.character(MITAB_polymer$edgetype), as.character(MITAB_polymer[,edge_weight_col])))
	} else {
		if (canonical_rep == "no") {
		edges_binary = unique(cbind(MITAB_binary$irogida, MITAB_binary$irogidb, as.character(MITAB_binary$experimental_role_A), as.character(MITAB_binary$experimental_role_B), MITAB_binary$irigid, as.character(MITAB_binary$edgetype), MITAB_binary[,edge_weight_col]))
		edges_complex = unique(cbind(MITAB_complex$irogida, MITAB_complex$irogidb, as.character(MITAB_complex$experimental_role_A), as.character(MITAB_complex$experimental_role_B), MITAB_complex$irigid, as.character(MITAB_complex$edgetype), MITAB_complex[,edge_weight_col]))
		edges_polymer = unique(cbind(MITAB_polymer$irogida, MITAB_polymer$irogidb, MITAB_polymer$irigid, as.character(MITAB_polymer$edgetype), as.character(MITAB_polymer[,edge_weight_col])))
		}
		if (canonical_rep == "yes") {
		edges_binary = unique(cbind(MITAB_binary$icrogida, MITAB_binary$icrogidb, as.character(MITAB_binary$experimental_role_A), as.character(MITAB_binary$experimental_role_B), MITAB_binary$icrigid, as.character(MITAB_binary$edgetype), MITAB_binary[,edge_weight_col]))
		edges_complex = unique(cbind(MITAB_complex$icrogida, MITAB_complex$icrogidb, as.character(MITAB_complex$experimental_role_A), as.character(MITAB_complex$experimental_role_B), MITAB_complex$icrigid, as.character(MITAB_complex$edgetype), MITAB_complex[,edge_weight_col]))
		edges_polymer = unique(cbind(MITAB_polymer$icrogida, MITAB_polymer$icrogidb, MITAB_polymer$icrigid, as.character(MITAB_polymer$edgetype), as.character(MITAB_polymer[,edge_weight_col])))
		}
	}

	# 3. Evaluate direction of binaries (note that there are no baits-preys for polymers):
	if (directionality=="directed" & dim(edges_binary)[1]>0) {
		new_edges_binary = NULL
		for (i in 1:dim(edges_binary)[1]) {
			if (edges_binary[i,3]=="MI:0496(bait)" & edges_binary[i,4]=="MI:0498(prey)") {
				new_edges_binary = rbind(new_edges_binary, edges_binary[i,])
			}
			if (edges_binary[i,3]=="MI:0498(prey)" & edges_binary[i,4]=="MI:0496(bait)") {
				tmp = c(edges_binary[i,2], edges_binary[i,1], edges_binary[i,3:7])
				new_edges_binary = rbind(new_edges_binary, tmp)
			}
			if (edges_binary[i,3]=="MI:0496(bait)" & edges_binary[i,4]=="MI:0496(bait)") {
				tmp1 = edges_binary[i,]
				tmp2 = c(edges_binary[i,2], edges_binary[i,1], edges_binary[i,3:7])
				new_edges_binary = rbind(new_edges_binary, tmp1)
				new_edges_binary = rbind(new_edges_binary, tmp2)
			}
		}
		print("Only rows with bait information are included in directed edgeLists")
	} else {
		new_edges_binary = edges_binary
	}

	if (dim(new_edges_binary)[1] == 1){
		edges_binary = cbind(new_edges_binary[1], new_edges_binary[2], new_edges_binary[5], new_edges_binary[6], new_edges_binary[7])
	} else {
		edges_binary = cbind(new_edges_binary[,1:2], new_edges_binary[,5:7])
	}

	# 4. Converting complexes to either matrix or spoke model:
	n = 0; edges_complex_list = list(); edges_complex_rep = NULL

	if (complex_rep == "spoke") {
		vector_of_complexes = unique(edges_complex[,1])
		for (i in vector_of_complexes) {
			vector_subunits = edges_complex[which(edges_complex[,1]==i), 2]
			vector_baits = edges_complex[which(edges_complex[,1]==i), 4]
			number_baits = length(vector_baits[which(vector_baits=="MI:0496(bait)")])
			irig = edges_complex[which(edges_complex[,1]==i), 5][1]
			edgeType = edges_complex[which(edges_complex[,1]==i), 6][1]
			weight = edges_complex[which(edges_complex[,1]==i), 7][1]
			if (number_baits == 1) {
				spoke_center = which(vector_baits=="MI:0496(bait)")
			}
			if (number_baits == 0) {
				spoke_center = 1	# First protein on list
			}
			if (number_baits > 1) {
				spoke_center = which(vector_baits=="MI:0496(bait)")[1]
			}
			if (length(vector_subunits) > 1) {
				spoke_rest = vector_subunits[-spoke_center]
				for (j in 1:length(spoke_rest)) {
					n = n + 1
					spoke_line = c(vector_subunits[spoke_center], spoke_rest[j], irig, edgeType, weight)
					edges_complex_list[[n]] = spoke_line
				}
			}
		}
		edges_complex_rep = do.call(rbind, edges_complex_list)
	}

	if (complex_rep == "matrix") {
		vector_of_complexes = unique(edges_complex[,1])
		for (i in 1:length(vector_of_complexes)) {
			vector_subunits = edges_complex[which(edges_complex[,1]==vector_of_complexes[i]), 2]
			irig = edges_complex[which(edges_complex[,1]==vector_of_complexes[i]), 5][1]
			edgeType = edges_complex[which(edges_complex[,1]==vector_of_complexes[i]), 6][1]
			weight = edges_complex[which(edges_complex[,1]==vector_of_complexes[i]), 7][1]
			if (length(vector_subunits) > 1) {
				tmp = t(combn(vector_subunits, 2))
				edges_complex_list[[i]] = cbind(tmp, rep(irig, dim(tmp)[1]), rep(edgeType, dim(tmp)[1]), rep(weight, dim(tmp)[1]))
			}
		}
		edges_complex_rep = do.call(rbind, edges_complex_list)
	}

	if (complex_rep == "bipartite" || dim(edges_complex)[1]==0) {
		edges_complex_rep = cbind(edges_complex[,1:2], edges_complex[,5:7])
	}

	# 5. Making unique list:
	total_unique = function(edgeList) {		# Remove repeated lines counting AB as equal to BA
		edgeList = unique(edgeList)
		ordered_edgeList = list()
		for (i in 1:dim(edgeList)[1]) {
			ordered_edgeList[[i]] = c(sort(edgeList[i,1:2]), edgeList[i,3:5])
		}
		edges = do.call(rbind, ordered_edgeList)
		edges = unique(edges)
	}

	if (multi_edge == "yes") {
		if (directionality=="directed") {
			final_list = unique(rbind(edges_binary, edges_complex_rep, edges_polymer))
		} else {
			final_list = total_unique(rbind(edges_binary, edges_complex_rep, edges_polymer))
		}
	} else {
		fused_list = NULL
		pre_final_list = rbind(edges_binary, edges_complex_rep, edges_polymer)
		tmp = data.frame(pre_final_list[,1], pre_final_list[,2])
		tmp2 = do.call(`paste`, c(unname(tmp), list(sep=".")))
		tmp3 = unique(tmp2)
		for (i in 1:length(tmp3)) {
			j = which(tmp2 == tmp3[i])
			tmp4 = cbind(pre_final_list[j,1], pre_final_list[j,2], rep(paste(sort(unique(pre_final_list[j,3])), sep="", collapse="|"), length(j)), rep(paste(sort(unique(pre_final_list[j,4])), sep="", collapse="|"), length(j)), pre_final_list[j,5])
			fused_list = rbind(fused_list, tmp4)
		}

		if (directionality=="directed") {
			final_list = unique(fused_list)
		} else {
			final_list = total_unique(fused_list)
		}
	}

	output = final_list
}
