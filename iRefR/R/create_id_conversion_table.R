####
# Generate a lookup table to convert any major protein ID to another protein ID:
####
create_id_conversion_table = function(MITAB_table, data_folder=getwd(), output_filename="id_conversion_table", IDs_to_include=c("uniprotkb", "refseq", "entrezgene/locuslink")) {

	# 1. Constructing output file names:
	if (data_folder == "data") {
		datafolder = system.file("data", package = "iRefR")
	} else if (data_folder == "home") {
		datafolder = R.home()
	} else {
		datafolder = data_folder
	}
	txtfile = paste(datafolder, "/", output_filename, ".txt", sep="")
	Rdatafile = paste(datafolder, "/", output_filename, ".RData", sep="")

	# 2. Generate table:
	cat("Generating ID Conversion table...\n")
	irogid = c(as.vector(MITAB_table$irogida), as.vector(MITAB_table$irogidb))
	icrogid = c(as.vector(MITAB_table$icrogida), as.vector(MITAB_table$icrogidb))
	alt = c(as.vector(MITAB_table$altA), as.vector(MITAB_table$altB))
	uid = c(as.vector(MITAB_table$X.uidA), as.vector(MITAB_table$uidB))
	initial_table = data.frame(irogid, icrogid, alt, uid)
	initial_table = unique(initial_table)
	extended_list = list(); n=0
	for (i in 1:dim(initial_table)[1]) {
		list_alts <- c(strsplit(as.character(initial_table[i,3]),"\\|")[[1]], as.character(initial_table[i,4]))
		for (j in 1:length(list_alts)) {
			tmp1 = strsplit(list_alts[j], "\\:")[[1]]
			if (setequal(IDs_to_include,"all")==TRUE) {
				if (tmp1[1]!="rogid" & tmp1[1]!="irogid") {
					n = n+1
					extended_list[[n]] = c(tmp1[1], tmp1[2], initial_table[i,1], initial_table[i,2])
				}
			} else {
				if (tmp1[1] %in% IDs_to_include) {
					n = n+1
					extended_list[[n]] = c(tmp1[1], tmp1[2], initial_table[i,1], initial_table[i,2])
				}
			}
		}
		if (as.integer(dim(initial_table)[1]/4) == i) {cat("25% completed...\n")}
		if (as.integer(dim(initial_table)[1]/2) == i) {cat("50% completed...\n")}
		if (as.integer(dim(initial_table)[1]*3/4) == i) {cat("75% completed...\n")}
	}
	# Convert list to matrix:
	extended_table = unique(do.call(rbind, extended_list))
	rownames(extended_table) = NULL; colnames(extended_table) = c("id_type", "id_value", "irogid", "icrogid")

	# 4. Generate output:
	write.table(extended_table, file = txtfile, sep="\t", row.names=FALSE, append=FALSE)
	cat("Conversion table has been copied to:\n")
	cat(paste(txtfile, "\n"))
	save(file = Rdatafile, list = "extended_table")

	output = extended_table
}
