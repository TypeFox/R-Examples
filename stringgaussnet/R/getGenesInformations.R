getGenesInformations <-
function(Identifiers,ensembl)
{
	if (requireNamespace("biomaRt",quietly=TRUE)) {Identifiers1=biomaRt::getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),filters="hgnc_symbol",values=Identifiers,mart=ensembl)} else {stop("biomaRt package must be installed to use this function")}
	Identifiers1=Identifiers1[which(!(duplicated(Identifiers1$hgnc_symbol))),]
	if (requireNamespace("biomaRt",quietly=TRUE)) {Identifiers2=biomaRt::getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),filters="ensembl_gene_id",values=Identifiers,mart=ensembl)} else {stop("biomaRt package must be installed to use this function")}
	Identifiers2=Identifiers2[which(!(duplicated(Identifiers2$ensembl_gene_id))),]
	TousIdentifiantsEnsembl=unique(c(Identifiers1$ensembl_gene_id,Identifiers2$ensembl_gene_id)) 
	
	if (requireNamespace("biomaRt",quietly=TRUE)) {GenesInformations=biomaRt::getBM(attributes=c("ensembl_gene_id","go_id"),filters="ensembl_gene_id",values=TousIdentifiantsEnsembl,mart=ensembl)} else {stop("biomaRt package must be installed to use this function")}
	if (requireNamespace("GO.db",quietly=TRUE)) {GoTermDesc<-GO.db::GOTERM} else {stop("GO.db package must be installed to use this function")}
	if (requireNamespace("AnnotationDbi",quietly=TRUE)) {Ontologies=AnnotationDbi::Ontology(GenesInformations$go_id)} else {stop("AnnotationDbi package must be installed to use this function")}
	if (requireNamespace("AnnotationDbi",quietly=TRUE)) {Terms=AnnotationDbi::Term(GenesInformations$go_id)} else {stop("AnnotationDbi package must be installed to use this function")}
	GOTerms=cbind(GenesInformations,Ontologies,Terms)
	GOTerms=GOTerms[which(GOTerms$Ontologies=="CC"),] 
	GOTerms=GOTerms[grep("(extracellular)|(^plasma membrane$)|(^cytosol$)|(^cytoplasm$)|(^nucleus$)|(^nucleoplasm$)",GOTerms$Terms),] 
	IdentifiantsLocalisations=list() 
	IdentifiantsLocalisations[["nuclear"]]=unique(GOTerms$ensembl_gene_id[grep("(^nucleus$)|(^nucleoplasm$)",GOTerms$Terms)]) 
	IdentifiantsLocalisations[["extracellular"]]=unique(GOTerms$ensembl_gene_id[grep("(extracellular)",GOTerms$Terms)]) 
	IdentifiantsLocalisations[["plasma membrane"]]=unique(GOTerms$ensembl_gene_id[grep("(^plasma membrane$)",GOTerms$Terms)]) 
	IdentifiantsLocalisations[["cytoplasm"]]=unique(GOTerms$ensembl_gene_id[grep("(^cytosol$)|(^cytoplasm$)",GOTerms$Terms)]) 
	GenesLocalisations=as.data.frame(cbind(TousIdentifiantsEnsembl,rep(NA,length(TousIdentifiantsEnsembl)))) 
	names(GenesLocalisations)=c("ensembl_gene_id","localization")
	for (Localisation in names(IdentifiantsLocalisations))
	{
		GenesLocalisations$localization[which(GenesLocalisations$ensembl_gene_id %in% IdentifiantsLocalisations[[Localisation]] & is.na(GenesLocalisations$localization))]=Localisation 
	}
	
	if (requireNamespace("biomaRt",quietly=TRUE)) {OtherInformations=biomaRt::getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "chromosome_name", "band", "strand", "start_position", "end_position", "description"),filters="ensembl_gene_id",values=TousIdentifiantsEnsembl,mart=ensembl)} else {stop("biomaRt package must be installed to use this function")}
	OtherInformations=OtherInformations[which(!duplicated(OtherInformations$hgnc_symbol)),] 
	AllOtherInformations=merge(GenesLocalisations,OtherInformations,by="ensembl_gene_id",all.x=T,sort=F) 
	
	rownames(AllOtherInformations)=AllOtherInformations$hgnc_symbol 
	rownames(AllOtherInformations)[which(AllOtherInformations$hgnc_symbol=="")]=AllOtherInformations$ensembl_gene_id[which(AllOtherInformations$hgnc_symbol=="")] 
	
	return(AllOtherInformations)
}
