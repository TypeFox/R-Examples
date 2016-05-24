# MAKING RNCBI PACKAGE
collapse <- function(...){paste(...,sep="",collapse="")}

EUtilsURL <- function(type="esearch",db="pubmed"){
	
	# CONSTRUCT ANY SERVICE TYPE AND DATABASE
	url <- 	"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/type.fcgi?db=DB&"
sub("(.*)(type)(.*)(DB)(.*)",collapse("\\1",type,"\\3",db,"\\5"),url)
}


EUtilsQuery <- function(query,type="esearch",db="pubmed",...){

	# CREATE URL WITH REQUEST
	PubMedURL <- EUtilsURL(type,db)
	Query <- gsub(" ","+",query)

	# PARTIAL MATCH OF POSSIBLE LIMITS
	OPTIONS <- c(
		"retstart",
		"retmax",
		"rettype",
		"field",
		"datetype",
		"reldate",
		"mindate",
		"maxdate"
	)
	
	ArgList <- list(...) # NEED FUNCTION PARTIAL MATCH HERE
	
	if(length(ArgList)==0){
		ArgList$retmax <- 1000 # DEFAULT 1000 RECORDS RETURNED
	}
	else{
		WhichArgs <- pmatch(names(ArgList),OPTIONS)	
		if(any(is.na(WhichArgs))||sapply(WhichArgs,length)>1)
			stop("Error in specified limits.")
		names(ArgList) <- OPTIONS[WhichArgs]
		if(all(names(ArgList)!="retmax"))
			ArgList$retmax <- 1000 # DEFAULT 1000 RECORDS RETURNED
	}
	
	ArgList$tool <- "RISmed"
	ArgList$email <- "s.a.kovalchik@gmail.com"
	
	# REPLACE RETMAX IF NOT USED
	ArgStr <- paste(names(ArgList),unlist(ArgList),sep="=")
	ArgStr <- paste(ArgStr,collapse="&")

paste(PubMedURL,"term=",Query,"&",ArgStr,sep="",collapse="")	
}

ParseTags <- function(lines){

	StripLines <- sapply(lines, function(x)gsub("> +",">",x)) # REPLACE WHITE SPACE
	
	Fields <- c("Count",
				"RetMax",
				"RetStart",
				"Id",
				"QueryTranslation")
			
	Patterns <- paste("(.*<",Fields,">)(.*)(<\\/",Fields,">).*",sep="")

	FieldIndex <- lapply(Patterns, function(pattern) grep(pattern, StripLines))
	FieldIndex[[1]] <- FieldIndex[[1]][1] # TAKE COUNT OF COMBINED QUERY

	Values <- lapply(1:length(Fields),
					function(i){
						result <- sapply(StripLines[FieldIndex[[i]]],
										function(x) sub(Patterns[i],"\\2",x),USE.NAMES=FALSE)
as.vector(result)
})


names(Values) <- Fields
Values$Count <- as.numeric(Values$Count)
Values$RetMax <- as.numeric(Values$RetMax)
Values$RetStart <- as.numeric(Values$RetStart)

Values
}


SplitIDs <- function(ids){

	if(length(ids)>200){
		group <- rep(1:ceiling(length(ids)/200),each=200)
		group <- group[1:length(ids)]
		split(ids,group)
	}
	else{
		list(ids)
	}

}


EUtilsGet <- function(x, type="efetch", db="pubmed"){

	if(class(x)[1]=="EUtilsSummary"){
		query <- x@querytranslation
		x <- x@PMID
	}
	else{
		query <- ""
	}
	
	IDList <- SplitIDs(x)
	Result <- lapply(IDList, EUtilsSubGet, type = type, db = db)
	Result <- unlist(Result, recursive=FALSE)

	if(type=="efetch"&db=="pubmed"){
		
		Result <- Medline(Result, query)
	}
			
Result
}


EUtilsSubGet <- function(ids, type="efetch", db="pubmed"){

	FetchURL <- EUtilsURL(type,db=db)
	IDStr <- collapse("id=",paste(ids,collapse=","))
	EUtilsFetch <- collapse(FetchURL,IDStr)	
	res <- readLines(collapse(EUtilsFetch,"&retmode=xml"), warn = FALSE, encoding = "UTF-8")
	
	pubstatus <- c(
		"received", 
		"accepted",
		"epublish",
		"ppublish",
		"pmc",
		"pubmed"
	)
	
	if(db=="pubmed"){
	
	ArticleList <- mapply(GroupArticle, start = ArticleStart(res),
                              end = ArticleEnd(res),
                              MoreArgs = list(.obj = res), 
                              SIMPLIFY = FALSE)

										
	ParseEUtilsFetch <- lapply(ArticleList, function(x){
		for(i in pubstatus){
			if(any(grepl("PubStatus", x) & grepl(i, x))){
				index <- which(grepl("PubStatus", x) & grepl(i, x))
				x[index + 1:5] <- gsub("(Year|Month|Day|Hour|Minute)",paste(c("\\1",First_Upper(i)),collapse = ""), x[index + 1:5])
			}
		}	

		if(any(grepl("AbstractText", x) & grepl("Label", x))){
			index <- which(grepl("AbstractText", x) & grepl("Label", x))
			labels <- sub("(.*Label.*[[:punct:]])([A-Z].*[A-Z])([[:punct:]].*Nlm.*)","\\2",x[index])
			for(i in 1:length(labels))
				x[index[i]] <- sub("(.*>)(.*>)",paste(c("\\1",labels[i],": \\2"), collapse = ""),x[index[i]])
		}
		
		val <- GetValues(x[LinesWithValues(x)])
		names(val) <- GetFields(x[LinesWithValues(x)])
	val
	})

	}
	else{
		tags <- GetFields(res[LinesWithValues(res)])
		ParseEUtilsFetch <- GetValues(res[LinesWithValues(res)])
		names(ParseEUtilsFetch) <- tags
	}
	
ParseEUtilsFetch
}

First_Upper <- function(x) paste(toupper(substr(x,1,1)),substr(x,2,nchar(x)),collapse = "", sep = "")

LinesWithValues <- function(.obj){
	grep(">(\\[?[a-zA-Z]|[0-9]).*<",.obj)
}

GetFields <- function(.obj){
	sub("(.*<)([a-zA-Z]+)(.*>)(\\[?([a-zA-Z]|[0-9]).*<..*>.*)","\\2",.obj)
}

GetValues <- function(.obj){
	sub("(.*<)([a-zA-Z]+)(.*>)(\\[?([a-zA-Z]|[0-9]).*)(<..*>.*)","\\4",.obj)
}

ArticleStart <- function(.obj) which(.obj=="<PubmedArticle>")
ArticleEnd <- function(.obj) which(.obj=="</PubmedArticle>")
GroupArticle <- function(start, end, .obj) .obj[start:end] 








