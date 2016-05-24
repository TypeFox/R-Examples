####
# Convert a table from complexList format to MITAB format:
####
convert_complexList_to_MITAB = function(complexList, MITAB_table, canonical_rep="yes", include_generated_complexes="no") {
	# 1. Separating known from generated complexes:
	found_dots = grep("\\.", complexList[,1])
	if (length(found_dots) > 0) {
		generated_complexes = complexList[found_dots,]
		known_complexes = complexList[-found_dots,]
	} else {
		known_complexes = complexList
		generated_complexes = NULL
	}

	# 2. Generating MITAB from known complexes:
	known_MITAB = NULL
	if (dim(known_complexes)[1] > 0) {
		list_rigids = known_complexes[,1]
		if (canonical_rep=="no") {
			known_MITAB = MITAB_table[which(MITAB_table$irigid %in% list_rigids),]
		} else {
			known_MITAB = MITAB_table[which(MITAB_table$icrigid %in% list_rigids),]
		}
	}

	# 3. Generating MITAB from generated complexes:
	MITAB_output = NULL
	if (include_generated_complexes=="yes") {
		MITAB_table = MITAB_table[which(MITAB_table$edgetype=="X"),]
		if(canonical_rep == "no") {
			table_a = data.frame(MITAB_table$sourcedb, MITAB_table$pmids, MITAB_table$method, MITAB_table$irogida)
			table_b = data.frame(MITAB_table$sourcedb, MITAB_table$pmids, MITAB_table$method, MITAB_table$irogidb)
		} else {
			table_a = data.frame(MITAB_table$sourcedb, MITAB_table$pmids, MITAB_table$method, MITAB_table$icrogida)
			table_b = data.frame(MITAB_table$sourcedb, MITAB_table$pmids, MITAB_table$method, MITAB_table$icrogidb)
		}
		MITAB_indexed_a = do.call(`paste`, c(unname(table_a), list(sep=".")))
		MITAB_indexed_b = do.call(`paste`, c(unname(table_b), list(sep=".")))
		tmp_a = MITAB_table[which(MITAB_indexed_a %in% generated_complexes),]
		tmp_b = MITAB_table[which(MITAB_indexed_b %in% generated_complexes),]
		MITAB_output = rbind(tmp_a, tmp_b)
	}

	output = unique(rbind(known_MITAB, MITAB_output))
}

		#MITAB_table = MITAB_table[which(MITAB_table$no baits
		#MITAB_list = list()
		#for (i in 1:dim(generated_complexes)[1]) {
			#tmp = strsplit(generated_complexes[i], "\\.")[[1]]
			#if (tmp[[4]] == "none") {
			#this_table = MITAB_table[which(MITAB_table$sourcedb==tmp[[1]] & MITAB_table$pmids==tmp[[2]] & MITAB_table$method==tmp[[3]]),]
			#MITAB_output = rbind(MITAB_output, this_table)
			#} else {
			#if (canonical_rep=="no") {
				#MITAB_list[[i]] = MITAB_table[which(MITAB_table$sourcedb==tmp[[1]] & MITAB_table$pmids==tmp[[2]] & MITAB_table$method==tmp[[3]] & ((MITAB_table$irogida==tmp[[4]] & MITAB_table$experimental_role_A=="MI:0496(bait)") | (MITAB_table$irogidb==tmp[[4]] & MITAB_table$experimental_role_B=="MI:0496(bait)"))),]
			#} else {
				#MITAB_list[[i]] = MITAB_table[which(MITAB_table$sourcedb==tmp[[1]] & MITAB_table$pmids==tmp[[2]] & MITAB_table$method==tmp[[3]] & ((MITAB_table$icrogida==tmp[[4]] & MITAB_table$experimental_role_A=="MI:0496(bait)") | (MITAB_table$icrogidb==tmp[[4]] & MITAB_table$experimental_role_B=="MI:0496(bait)"))),]
			#}
			#}
		#}
		#MITAB_output = do.call(rbind, MITAB_list)
