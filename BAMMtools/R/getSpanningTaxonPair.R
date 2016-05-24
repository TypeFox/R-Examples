#############################################################
#
#	getSpanningTaxonPair(....)
#
#	returns pair of taxa that span a given taxon set

getSpanningTaxonPair <- function(phy, taxset){
	
	if (! sum(taxset %in% phy$tip.label) > 0){
		cat('Some species in taxset that are not in tree\n');
		taxset <- taxset[taxset %in% phy$tip.label];
	}
	
	dt <- drop.tip(phy, setdiff(phy$tip.label, taxset));
	
	return(c(dt$tip.label[1], dt$tip.label[length(dt$tip.label)]));
}
