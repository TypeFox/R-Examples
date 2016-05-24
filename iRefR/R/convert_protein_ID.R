####
# Convert protein IDs between iROGs, icROGs, geneIDs, RefSeqs and UniProtIDs:
####
convert_protein_ID = function(input_id_type, input_id_value, output_id_type, id_conversion_table) {

	# 1. Get icrogid (only one per input value):
	if (input_id_type == "irogid") {
		protein_icrogid = id_conversion_table[which(id_conversion_table[,3]==input_id_value), 4][1]
	}
	if (input_id_type == "icrogid" ) {protein_icrogid = input_id_value}
	if (input_id_type == "uniprotkb" || input_id_type == "refseq" || input_id_type == "entrezgene/locuslink") {
		protein_icrogid = id_conversion_table[which(id_conversion_table[,1]==input_id_type & id_conversion_table[,2]==input_id_value), 4]
	}

	# 2. Find output (one or more per icrogid):
	if (output_id_type == "irogid") {
		output_id_value = id_conversion_table[which(id_conversion_table[,4]==protein_icrogid), 3]
	}
	if (output_id_type == "icrogid") {output_id_value = protein_icrogid}
	if (output_id_type == "uniprotkb" || output_id_type == "refseq" || output_id_type == "entrezgene/locuslink") {
		output_id_value = id_conversion_table[which(id_conversion_table[,1]==output_id_type & id_conversion_table[,4]==protein_icrogid), 2]
	}

	output_id_value = unique(output_id_value)
	output_id_value
}
