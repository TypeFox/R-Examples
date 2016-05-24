####
# Get statistical information about some of the columns of a given iRefIndex/MITAB table:
####
summary_table = function(MITAB_table) {

	# Counts:
	number_irigids = length(unique(MITAB_table$irigid))
	number_icrigids = length(unique(MITAB_table$icrigid))
	number_irogids = length(unique(c(MITAB_table$irogida, MITAB_table$irogidb)))
	number_icrogids = length(unique(c(MITAB_table$icrogida, MITAB_table$icrogidb)))
	number_pmids = length(unique(MITAB_table$pmids))
	number_methods = length(unique(MITAB_table$method))

	# Percentages:
	db_table = table(MITAB_table$sourcedb)
	int_type_table = table(MITAB_table$edgetype)
	particip_table = table(MITAB_table$numParticipants)

	# Output:
	output = list(number_non_canonical_interactions=number_irigids, number_canonical_interactions=number_icrigids, number_non_canonical_proteins=number_irogids, number_canonical_proteins=number_icrogids, number_publications=number_pmids, number_experimental_methods=number_methods, source_db_distribution=db_table, interaction_type_distribution=int_type_table, number_participants_distribution=particip_table)

}
