####
# Get information on interactors and complexes, for a given protein:
####
summary_protein = function(id_type, id_value, MITAB_table) {
	# 1. iROG ID:
	if (id_type == "irogid") {
		all_interactions_A = MITAB_table[which(MITAB_table$irogida == id_value),]
		all_interactions_B = MITAB_table[which(MITAB_table$irogidb == id_value),]
		all_interactions = unique(rbind(all_interactions_A, all_interactions_B))

		all_degree = length(unique(as.vector(all_interactions$irigid)))
		interactions_X = all_interactions[which(all_interactions$edgetype == "X"),]
		number_interactions_X = length(unique(as.vector(interactions_X$irigid)))
		interactors_X = sort(unique(c(as.vector(interactions_X$irogida), as.vector(interactions_X$irogidb))))
		interactors_X =  interactors_X[which(interactors_X != "-" & interactors_X != id_value)]
		
		interaction_this_value_C = all_interactions[which(all_interactions$edgetype == "C"),]
		irigid_this_C = unique(as.vector(interaction_this_value_C$irigid))
		vector_interactors_C = NULL
		for (i in irigid_this_C) {
			interactions_C = MITAB_table[which(MITAB_table$irigid == i),]
			interactors_C = sort(unique(as.vector(interactions_C$irogidb)))
			interactors_C =  interactors_C[which(interactors_C != "-" & interactors_C != id_value)]
			interactors_C = paste(interactors_C, sep="", collapse=",")
			vector_interactors_C = c(vector_interactors_C, interactors_C)
		}
		number_interactions_C = length(irigid_this_C)
		interactors_C = vector_interactors_C

		interactions_Y = all_interactions[which(all_interactions$edgetype == "Y"),]
		number_interactions_Y = length(unique(as.vector(interactions_Y$irigid)))
		interactors_Y = sort(unique(c(as.vector(interactions_Y$irogida), as.vector(interactions_Y$irogidb))))
		interactors_Y =  interactors_Y[which(interactors_Y != "-" & interactors_Y != id_value)]
	}

	# 2. icROG ID:
	if (id_type == "icrogid") {
		all_interactions_A = MITAB_table[which(MITAB_table$icrogida == id_value),]
		all_interactions_B = MITAB_table[which(MITAB_table$icrogidb == id_value),]
		all_interactions = unique(rbind(all_interactions_A, all_interactions_B))

		all_degree = length(unique(as.vector(all_interactions$icrigid)))
		interactions_X = all_interactions[which(all_interactions$edgetype == "X"),]
		number_interactions_X = length(unique(as.vector(interactions_X$icrigid)))
		interactors_X = sort(unique(c(as.vector(interactions_X$icrogida), as.vector(interactions_X$icrogidb))))
		interactors_X =  interactors_X[which(interactors_X != "-" & interactors_X != id_value)]
		
		interaction_this_value_C = all_interactions[which(all_interactions$edgetype == "C"),]
		irigid_this_C = unique(as.vector(interaction_this_value_C$icrigid))
		vector_interactors_C = NULL
		for (i in irigid_this_C) {
			interactions_C = MITAB_table[which(MITAB_table$icrigid == i),]
			interactors_C = sort(unique(as.vector(interactions_C$icrogidb)))
			interactors_C =  interactors_C[which(interactors_C != "-" & interactors_C != id_value)]
			interactors_C = paste(interactors_C, sep="", collapse=",")
			vector_interactors_C = c(vector_interactors_C, interactors_C)
		}
		number_interactions_C = length(irigid_this_C)
		interactors_C = vector_interactors_C

		interactions_Y = all_interactions[which(all_interactions$edgetype == "Y"),]
		number_interactions_Y = length(unique(as.vector(interactions_Y$icrigid)))
		interactors_Y = sort(unique(c(as.vector(interactions_Y$icrogida), as.vector(interactions_Y$icrogidb))))
		interactors_Y =  interactors_Y[which(interactors_Y != "-" & interactors_Y != id_value)]
	}

	# 3. Generating output:
	if (dim(all_interactions)[1] == 0) {
		table_protein_info = "ERROR. Protein ID not found."
		interactors_X =  "ERROR. Protein ID not found."
		interactors_C =  "ERROR. Protein ID not found."
		interactors_Y =  "ERROR. Protein ID not found."
	} else {
		if (id_type == "irogid") {
			table_names = c(id_type, "number_interactions_X(#irigids)", "number_interactions_C(#irigids)", "number_interactions_Y(#irigids)", "total_degree(#irigs)")
		} else {
			table_names = c(id_type, "number_interactions_X(#icrigids)", "number_interactions_C(#icrigids)", "number_interactions_Y(#icrigids)", "total_degree(#icrigs)")
		}
		table_values = c(id_value, number_interactions_X, number_interactions_C, number_interactions_Y, all_degree)
		table_protein_info = data.frame(table_names, table_values)
	}

	output = list(summary_interactions_per_protein=table_protein_info, interactors_binary=interactors_X, interactors_complex=interactors_C, interactors_polymer=interactors_Y)

}
