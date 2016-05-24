####
# Convert a table from MITAB format to complexList format:
####
convert_MITAB_to_complexList = function(MITAB_table, canonical_rep="yes", include_generated_complexes="no", list_methods="default", bait_use="only_baits", node_names="rogs") {

	# 1. Separate complexes from binaries:
	MITAB_binary = select_interaction_type("binary", MITAB_table)
	MITAB_complex = select_interaction_type("complex", MITAB_table)

	# 2. Generating complexList from known complexes:
	complex_name = list(); complexList = list()
	if (node_names == "uids") {
		irigids = unique(MITAB_complex$irigid)
		for (i in irigids) {
			complex_name[[i]] = i
			this_complex = sort(unique(as.character(MITAB_complex[which(MITAB_complex$irigid==i), "uidB"])))
			complexList[[i]] = paste(this_complex, collapse=",", sep="")
		}
	} else {
		if (canonical_rep=="no") {
			irigids = unique(MITAB_complex$irigid)
			for (i in irigids) {
				complex_name[[i]] = i
				#complex_name[[i]] = unique(MITAB_complex[which(MITAB_complex$irigid==i), "irogida"])[1]	# this [1] for cases where one icrigid = more than one icrogida, which is strange
				this_complex = sort(unique(MITAB_complex[which(MITAB_complex$irigid==i), "irogidb"]))
				complexList[[i]] = paste(this_complex, collapse=",", sep="")
			}
		}
		if (canonical_rep == "yes") {
			icrigids = unique(MITAB_complex$icrigid)
			for (i in icrigids) {
				complex_name[[i]] = i
				this_complex = sort(unique(MITAB_complex[which(MITAB_complex$icrigid==i), "icrogidb"]))
				complexList[[i]] = paste(this_complex, collapse=",", sep="")
			}
		}
	}
	complexName_table = do.call(rbind, complex_name)
	complexList_table = do.call(rbind, complexList)
	complexList_table = data.frame(complexName_table, complexList_table)

	# 3. Generating possibly binary-represented complexes from binaries:
	generate_complexes_from_binary_represented = function(MITAB_binary, canonical_rep, list_methods, node_names) {
		# Generating complexes from AP/coip-like data reported as binaries:
		# Here, a complex is not something that copurifies and acts as a functional unit; here,
		# everything that has 3 or more members, copurifies through the same method, in the same
		# pmid, in the same DB, and has information about at least one bait and at least two
		# preys, is assumed to be a binary-represented complex.

		# 3.1. Creating a table of records with copurification-method information:
		# List of selected AP/coip-like methods:
		if (list_methods=="default") {
			list_of_methods = c("MI:0004\\(", "MI:0006\\(", "MI:0007\\(", "MI:0019\\(", "MI:0027\\(", "MI:0028\\(", "MI:0029\\(", "MI:0059\\(", "MI:0061\\(", "MI:0071\\(", "MI:0096\\(pull down\\)", "MI:0114\\(", "MI:0226\\(", "MI:0227\\(", "MI:0401\\(", "MI:0437\\(", "MI:0676\\(", "MI:0858\\(", "MI:0963\\(")
		} else {
			list_of_methods = paste(list_methods, "\\(", sep="")
		}

		binaries_with_method = NULL
		for (i in list_of_methods) {
			index_misplaced = grep(i, MITAB_binary$method)
			binaries_with_method = rbind(binaries_with_method, MITAB_binary[index_misplaced,])
		}

		# 3.2. Generating complexes with the above-specified conditions:
		list_groups_method_pmid_db = unique(data.frame(binaries_with_method$method, binaries_with_method$pmids, binaries_with_method$sourcedb))

		list_complexes = list(); list_complexes_names = list(); n = 0
		for (i in 1:dim(list_groups_method_pmid_db)[1]) {
			tmp = binaries_with_method[which(binaries_with_method$method==list_groups_method_pmid_db[i,1] & binaries_with_method$pmids==list_groups_method_pmid_db[i,2] & binaries_with_method$sourcedb==list_groups_method_pmid_db[i,3]),]
			if (node_names == "uids") {
				list_prots = c(tmp$X.uidA, tmp$uidB)
			} else {
				if (canonical_rep=="no") {
					list_prots = c(tmp$irogida, tmp$irogidb)
				} else {
					list_prots = c(tmp$icrogida, tmp$icrogidb)
				}
			}
			list_roles = c(as.character(tmp$experimental_role_A), as.character(tmp$experimental_role_B))
			baits = unique(list_prots[which(list_roles=="MI:0496(bait)")])
			if (length(baits)>0) {
				for (j in baits) {
					if (node_names == "uids") {
						preys = sort(unique(c(tmp[which(tmp$X.uidA==j),"uidB"], tmp[which(tmp$uidB==j),"X.uidA"])))
					} else {
						if (canonical_rep=="no") {
							preys = sort(unique(c(tmp[which(tmp$irogida==j),"irogidb"], tmp[which(tmp$irogidb==j),"irogida"])))
						} else {
							preys = sort(unique(c(tmp[which(tmp$icrogida==j),"icrogidb"], tmp[which(tmp$icrogidb==j),"icrogida"])))
						}
					}
					if (length(preys)>1) {
						n = n + 1
						list_complexes[[n]] = paste(c(j, preys), collapse=",", sep="")
						list_complexes_names[[n]] = paste(unique(tmp$sourcedb), unique(tmp$pmids), unique(tmp$method), j, sep=".")
					}
				}
			} else {
				if (bait_use=="include_non_baits" & length(unique(list_prots))>3) {
					candidate_centers = unique(list_prots)
					for (j in candidate_centers) {
						if (node_names == "uids") {
							tmp2 = c(tmp[which(tmp$X.uidA==j),"uidB"], tmp[which(tmp$uidB==j),"X.uidA"])
						} else {
							if (canonical_rep=="no") {
								tmp2 = c(tmp[which(tmp$irogida==j),"irogidb"], tmp[which(tmp$irogidb==j),"irogida"])
							} else {
								tmp2 = c(tmp[which(tmp$icrogida==j),"icrogidb"], tmp[which(tmp$icrogidb==j),"icrogida"])
							}
						}
						tmp3 = sort(unique(c(j, tmp2)))
						if (length(tmp3)>3) {
							n = n + 1
							list_complexes[[n]] = paste(tmp3, collapse=",", sep="")
							list_complexes_names[[n]] = paste(unique(tmp$sourcedb), unique(tmp$pmids), unique(tmp$method), j, sep=".")
						}
					}
				}
			}
		}

		column_ids = do.call(rbind, list_complexes_names)
		column_subunits = do.call(rbind, list_complexes)
		the_table = cbind(column_ids, column_subunits)

		rownames(the_table) = NULL
		if (node_names == "uids") {
			colnames(the_table) = c("complex_ID", "UID_subunits")
		} else {
			if (canonical_rep=="no") {
				colnames(the_table) = c("complex_ID", "iROG_ID_subunits")
			} else {
				colnames(the_table) = c("complex_ID", "icROG_ID_subunits")
			}
		}

		output = the_table
	}

	# 3.3. Using the previous function to generate complexes:
	if (include_generated_complexes=="yes") {
		generated_table = generate_complexes_from_binary_represented(MITAB_binary, canonical_rep, list_methods, node_names)
		complexList_table = rbind(complexList_table, generated_table)
	}

	# 4. Output:
	output = complexList_table
}
