cited_function <- function(object){

if(class(object)[1] == "EUtilsSummary" && object@db != "pubmed")
	stop("Cited is only available for queries of the Pubmed database")

f <- function(id){

	base <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?retmode=xml&dbfrom=pubmed&id=ID&cmd=neighbor"
	base <- sub("ID", id, base)
	
	the_url <- url(description = base)
	on.exit(close(the_url))

	lines <- readLines(the_url)
	
	citedstart <- grep("pubmed_pubmed_citedin", lines)

	if(length(citedstart) == 0){
		return(0)
		}
	else{
		citedend <- grep("pubmed_pubmed", lines)
		citedend <- min(citedend[citedend > citedstart])
	
		tags <- grep("<Id>", lines)
		tags <- tags[tags > citedstart & tags < citedend]
	
		if(any(tags)){
			hits <- sub("(.*<ID>)([0-9]+)(.*)","\\2",lines[tags])
			hits <- unique(hits)
		length(hits[hits != id])
		}
		else{
			0
		}	
	}
}

sapply(FUN = f, object@PMID)
}
