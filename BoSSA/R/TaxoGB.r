`TaxoGB` <-
function(X,tsleep=3,organism="viridiplantae")
{
	query <- gsub(" ","+",paste(X,"+AND+",organism,"[organism]",sep=""))
	GBFile <- scan(paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?&db=taxonomy&retmax=20&term=",query,sep=""),what="",sep="\n",quiet=TRUE)
	GB <- xml(c("<root>",GBFile[3:length(GBFile)],"</root>"))
	ID <- NULL
	Family <- NULL
	Genus <- NULL
	subgenus <- NULL
	error <- GB['eSearchResult/ErrorList']
	if(length(error)==0)
	{
	    ID <- GB['eSearchResult/IdList/Id/#'][1]
	    GBFile2 <- scan(paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=",ID,"&retmode=xml",sep=""),what="",sep="\n",quiet=TRUE)
	    GB2 <- xml(c("<root>",GBFile2[3:length(GBFile2)],"</root>"))
	    Family <- GB2['TaxaSet/Taxon/LineageEx/Taxon/ScientificName/#'][GB2['TaxaSet/Taxon/LineageEx/Taxon/Rank/#']=="family"]
	    Genus <- GB2['TaxaSet/Taxon/LineageEx/Taxon/ScientificName/#'][GB2['TaxaSet/Taxon/LineageEx/Taxon/Rank/#']=="genus"]
	    subgenus <- GB2['TaxaSet/Taxon/LineageEx/Taxon/ScientificName/#'][GB2['TaxaSet/Taxon/LineageEx/Taxon/Rank/#']=="subgenus"]
	}
	CorrectedQ <- NULL
	if(length(error)>0)
	{
	  GBFile3 <- scan(paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/espell.fcgi?&db=pubmed&retmode=xml&term=",query,sep=""),what="",sep="\n",quiet=TRUE)
	  GB3 <- xml(c("<root>",GBFile3[3:length(GBFile3)],"</root>"))
	  CorrectedQ <- gsub("  and viridiplantae organism","",GB3['eSpellResult/CorrectedQuery/#'])
	}
	Sys.sleep(tsleep)
	out <- paste(X,ID,Family,Genus,subgenus,CorrectedQ,sep="\t")
	out
}

