concordanceFun <-
function(nuc, aa, nuclength = 648, aalength = 216, aminoAcid_Modal){
	VLF_match <- find.matching(nuc, aa, nuclength, aalength)

	position_matchingNuc <- nucleotide.matching.positions(VLF_match[[2]], nuclength)

	position_matchingAA <- aminoAcid.matching.positions(VLF_match[[1]], aalength)

	matching_comparison <- overall.matched(position_matchingNuc, position_matchingAA, nuclength, aalength)

	matching_codons <- matched.codon.position(matching_comparison)

	concordant_aaType_change <- concordant.to.modalchanges(matching_comparison, aminoAcid_Modal)

	All_aaType_change <- aaVLFs.to.modalchanges(aminoAcid_Modal, aa, aalength)
	
	foo <- list(matched = matching_comparison, codons = matching_codons, concordantType = concordant_aaType_change, aminoAcidType = All_aaType_change, concordNuc = length(matching_comparison), concordAA = sum(concordant_aaType_change), sequences = nrow(VLF_match[[1]]))
}
