setClass("EUtilsSummary",
	representation(
	    db = "character",
		count = "numeric",
		retmax = "numeric",
		retstart = "numeric",
		PMID = "character",
		querytranslation = "character"
	)
)

EUtilsSummary <- function(query,type="esearch",db="pubmed",url=NULL,encoding="unknown",...){

	if(is.null(url)){
		url <- EUtilsQuery(query,type,db,...)
	}

	lines <- readLines(url,warn=FALSE,encoding=encoding)
	res <- ParseTags(lines)
	
	EMPTYCHECK <- length(grep("<eSearchResult><Count>0<\\/Count>", lines))!=0
	
	if(EMPTYCHECK){
			res$Id <- character(0)
			res$Count <- 0
		}
	
	new("EUtilsSummary",
		db = db,
		count = res$Count,
		retstart = res$RetStart,
		retmax = res$RetMax,
		PMID = res$Id,
		querytranslation = res$QueryTranslation
		)
}

setMethod("print","EUtilsSummary",function(x,...) print(x@querytranslation))
setMethod("show","EUtilsSummary",function(object) print(object@querytranslation))
setMethod("summary","EUtilsSummary",function(object,...){

			cat("Query:\n")
			cat(object@querytranslation,"\n\n")
			cat("Result count: ",object@count)
			
			invisible(object@PMID)
			
})

# GENERICS
setMethod("QueryCount","EUtilsSummary",function(object) object@count)
setMethod("QueryId","EUtilsSummary",function(object) object@PMID)
setMethod("QueryTranslation","EUtilsSummary",function(object) object@querytranslation)
setMethod("Cited", "EUtilsSummary", cited_function)








