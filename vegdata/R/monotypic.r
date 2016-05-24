monotypic <- function(refl, nr.member = 1, reflist.type = c('Turboveg', 'EDIT'), write = FALSE, filename,  tv_home, ...) {
 reflist.type <- match.arg(reflist.type)	
 if(reflist.type=='Turboveg' & missing(tv_home)) tv_home <- tv.home()
 taxa <- load.taxlist(refl = refl, reflist.type=reflist.type , detailed = TRUE, ...)
 names(taxa) <- TCS.replace(names(taxa))
 AG <- table(taxa$IsChildTaxonOfID)
 AG <- names(AG[AG == nr.member])
 if(reflist.type == 'Turboveg') AG <- as.numeric(AG)
# names(taxa)[names(taxa)=='TAXONRANK'] <- 'TaxonRank'
 mono <- data.frame(Parent_NR = AG, Parent_Name=taxa$TaxonName[match(AG, taxa$TaxonUsageID)], Parent_Rank = taxa$TaxonRank[match(AG, taxa$TaxonUsageID)], MEMBER_NR=taxa$TaxonUsageID[match(AG,taxa$IsChildTaxonOfID)], MEMB_NAME=taxa$TaxonName[match(AG,taxa$IsChildTaxonOfID)], MEMB_Rank=taxa$TaxonRank[match(AG,taxa$IsChildTaxonOfID)])
 NbChildren <- as.integer(sapply(mono$MEMBER_NR, function(x) nrow(child(x, refl=refl, reflist.type=reflist.type, quiet=TRUE)) ), use.names=FALSE)
mono$NbChildren <- NbChildren
# AGG_NR,N,8,0  AGG_taxonR,C,13	MEMBER_NR,N,9,0	MEMB_NAME,C,67	MEMB_taxon,C,14
# mono$NbChildren <- unlist(mono$NbChildren)
#mono$NbChildren ch <- child(x, refl=refl, reflist.type=reflist.type) 
# if(write) write.csv2(mono[NbChildren == nr.member,], filename, row.names = FALSE)
 if(write) write.csv2(mono, filename, row.names = FALSE)
invisible(mono)
}

tv.mono <- function(...) stop('Deprecated. Use monotypic() instead.')
