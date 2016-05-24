extract_carriers = function(ped,site,fam,type="alleles",minor.allele=2)
# ped : pedigree coded in a ped file with either two alleles per variant ("alleles"), or a count of one allele ("count")
# site : site where to record carriers
# fam : ID of the family for which to extract carriers
{
if (!(type %in% c("alleles","count"))) stop ("Invalid type ",type)

if (type == "alleles")	
{	
if (missing(fam))
{
	if (length(unique(ped[,1]))>1) stop ("More than one family in ped data and no family specified.")
	genodat = ped[,5:6+2*site]
	if (any(genodat==0)) stop ("Alleles can't be labelled 0.")
	subid = ped[,2]
}	
else
{
	genodat = ped[ped[,1]==fam&ped[,6]==2,5:6+2*site]
	subid = ped[ped[,1]==fam&ped[,6]==2,2]		
}

carriers.bool = apply(genodat,1,function(vec) ifelse(any(is.na(vec)),FALSE,any(vec==minor.allele)) )
}
else # type == "count"
{
if (missing(fam))
{
	if (length(unique(ped[,1]))>1) stop ("More than one family in ped data and no family specified.")
	allelecount = ped[ped[,6]==2,6+site]
	subid = ped[,2]
}
else
{
	allelecount = ped[ped[,1]==fam&ped[,6]==2,6+site]
	subid = ped[ped[,1]==fam&ped[,6]==2,2]		
}
    carriers.bool = ifelse(is.na(allelecount),FALSE,allelecount>0)
}
# Return list of carriers
as.character(subid[carriers.bool])
}

