`InfoGB` <-
function(X,tsleep=3)
{
	GBID <- X
	GBFile <- scan(paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",GBID,"&rettype=gbc",sep=""),what="",sep="\n",quiet=TRUE)
	GB <- xml(c("<root>",GBFile[3:length(GBFile)],"</root>"))
	
	DateSub <- GB['INSDSet/INSDSeq/~^INSDSeq_create/#']
	Organism <- GB['INSDSet/INSDSeq/INSDSeq_organism/#']
	Taxo <- GB['INSDSet/INSDSeq/INSDSeq_taxonomy/#']
	nbauthor <- length(GB['INSDSet/INSDSeq/INSDSeq_references/INSDReference/INSDReference_authors'][[1]])-2
	authors <- NULL
	for (i in 1:nbauthor)
	{
	authors <- c(authors,GB[paste('INSDSet/INSDSeq/INSDSeq_references/INSDReference/INSDReference_authors/INSDAuthor[',i,']/#',sep="")][[1]])
	}
	authors <- paste(authors,collapse=" ; ")
	Titre <- GB['INSDSet/INSDSeq/INSDSeq_references/INSDReference/INSDReference_title/#'][[1]]
	Journal <- GB['INSDSet/INSDSeq/INSDSeq_references/INSDReference/INSDReference_journal/#'][[1]]
	pubmed <- GB['INSDSet/INSDSeq/INSDSeq_references/INSDReference/INSDReference_pubmed/#'][[1]]
	pubmedURL <- ""
	if(length(pubmed)>0) pubmedURL <- paste("http://www.ncbi.nlm.nih.gov/pubmed/",pubmed,sep="")
	isolate <- GB['INSDSet/INSDSeq/~^INSDSeq_feature/INSDFeature/INSDFeature_quals/INSDQualifier/INSDQualifier_value/#'][GB['INSDSet/INSDSeq/~^INSDSeq_feature/INSDFeature/INSDFeature_quals/INSDQualifier/INSDQualifier_name/#']=="isolate"]
	host <- GB['INSDSet/INSDSeq/~^INSDSeq_feature/INSDFeature/INSDFeature_quals/INSDQualifier/INSDQualifier_value/#'][GB['INSDSet/INSDSeq/~^INSDSeq_feature/INSDFeature/INSDFeature_quals/INSDQualifier/INSDQualifier_name/#']=="host"]
	Location <- GB['INSDSet/INSDSeq/~^INSDSeq_feature/INSDFeature/INSDFeature_quals/INSDQualifier/INSDQualifier_value/#'][GB['INSDSet/INSDSeq/~^INSDSeq_feature/INSDFeature/INSDFeature_quals/INSDQualifier/INSDQualifier_name/#']=="country"]
	DateEch <- GB['INSDSet/INSDSeq/~^INSDSeq_feature/INSDFeature/INSDFeature_quals/INSDQualifier/INSDQualifier_value/#'][GB['INSDSet/INSDSeq/~^INSDSeq_feature/INSDFeature/INSDFeature_quals/INSDQualifier/INSDQualifier_name/#']=="collection_date"]
	GPS <- GB['INSDSet/INSDSeq/~^INSDSeq_feature/INSDFeature/INSDFeature_quals/INSDQualifier/INSDQualifier_value/#'][GB['INSDSet/INSDSeq/~^INSDSeq_feature/INSDFeature/INSDFeature_quals/INSDQualifier/INSDQualifier_name/#']=="lat_lon"]
	source <- GB['INSDSet/INSDSeq/~^INSDSeq_feature/INSDFeature/INSDFeature_quals/INSDQualifier/INSDQualifier_value/#'][GB['INSDSet/INSDSeq/~^INSDSeq_feature/INSDFeature/INSDFeature_quals/INSDQualifier/INSDQualifier_name/#']=="isolation_source"]
	if(length(host)>0)
	{
		hostname <- host
		taxoHost <- TaxoGB(hostname)
	}
	else if(length(host)>0)
	{		
		hostname <- source
		taxoHost <- TaxoGB(hostname)
	}
	else taxoHost <- "\t\t\t\t\t"

	Sys.sleep(tsleep)

	out <- paste(X,Organism,isolate,Taxo,DateSub,DateEch,taxoHost,Location,GPS,authors,Titre,Journal,pubmedURL,sep="\t")
	out
}




