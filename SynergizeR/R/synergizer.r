#' Translate a set of biological identifiers into an selected alternative.
#'
#' This function will translate between sets of biological identifiers.
#'
#' @param authority A character containing any authoritative sources of identifier-mapping information.
#' @param species A character containing the Species. Note that the range of species supported 
#' depends on the choice of authority. Examples: Homo sapiens, Mus musculus.
#' @param domain This is the "namespace" (naming scheme) of the database identifiers the user wishes to translate.
#' Examples: embl, ipi
#' @param range This is the "namespace" (naming scheme) to which the user wishes to translate the input identifiers.
#' Examples: embl, ipi
#' @param ids a vector containing the ids to be translated
#' @param file NULL or a string containing the name of the file where the ids will be saved
#'
#' @return A vector containing the translated ids.
#'
#' @references http://llama.mshri.on.ca/synergizer/translate/
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' library('SynergizerR')
#' symbols.ids <- synergizer( authority = "ensembl", species = "Homo sapiens", domain="affy_hg_u95av2", range="hgnc_symbol",ids=c("1939_at","1503_at","1454_at") )
#'
#' entrez.ids <- synergizer( authority = "ensembl", species = "Homo sapiens", domain="hgnc_symbol", range="entrezgene",ids=c("snph", "pja1", "prkdc", "RAD21L1", "Rorc", "kcnk16") ) }
#'
synergizer <- function( authority = "ensembl", 
						species = "Homo sapiens", 
						domain = "hgnc_symbol", 
						range = "entrezgene", 
						ids = NULL, 
						file = NULL)
{
	if(is.null(ids)) stop("Insert a character vector with the ids to be translated!")
	else {}
	Sys.sleep(3) # In order to not have user host machine banned from using the service.
	t <- basicTextGatherer()
	h <- basicHeaderGatherer()

	curlPerform( url = "http://llama.mshri.on.ca/cgi/synergizer/serv",
				 .opts = list(httpheader = c('Content-Type' = "application/json"), verbose = FALSE),
				 curl = getCurlHandle(),
				 postfields = toJSON( list( method="translate",params=list(list(authority=authority, species=species, 
				 domain=domain, range=range, ids=ids)), id=123 ) ),
				 writefunction = t$update,
				 headerfunction = h$update
				) 

	output <- list(data = t$value(),
	               status = h$value()[['status']],
	               status.message = h$value()[['statusMessage']])
	# check http status
	httpstatus <- as.numeric(output$status)
	if (httpstatus != 200) {
		print(output$status.message); break
		}
	else response <- fromJSON(output$data, nullValue = NA) # map NULL JSON 'null' value to R 'NA'
	# return(response$result) # it outputs the original list
	lens <- sapply(response$result, length)
	output <- matrix(nrow=length(response$result), ncol=2)
	for (i in 1:length(response$result)) {
		output[i,1] = response$result[[i]][1]
		output[i,2] = paste(response$result[[i]][2:lens[i]],collapse="|")
	}
	colnames(output) <- c( deparse(substitute(domain)), deparse(substitute(range)) ) 
	if(!is.null(file)) write.table( output, file = file, quote = F, sep = "\t", row.names = FALSE, col.names = TRUE )
	# return( list(as.data.frame(output, stringsAsFactors=FALSE), response$result) ) # debug
	return( as.data.frame(output, stringsAsFactors=FALSE) )
}
