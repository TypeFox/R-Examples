#' Returns current version of the Synergizer server as a string.
#' 
#' @return A one-element character vector containing current version of the Synergizer server
#'
#' @references http://llama.mshri.on.ca/synergizer/translate/
#'
#' @export
#'
#' @examples
#' \dontrun{ 
#' library('SynergizerR')
#' server_version() 
#' }
#'
server_version <- function(){
	Sys.sleep(3) # In order to not have user host machine banned from using the service.
	header.field = c('Content-Type' = "application/json")
	curl <- getCurlHandle() 
	curlSetOpt(.opts = list(httpheader = header.field, verbose = FALSE), curl = curl) 
	t <- basicTextGatherer()
	h <- basicHeaderGatherer()

	body = curlPerform(url = "http://llama.mshri.on.ca/cgi/synergizer/serv",
	                   curl = curl,
	                   postfields = '{"method":"version","params":[],"id":0}',
	                   writefunction = t$update,
	                   headerfunction = h$update)	
	
	output <- list(data = t$value(),
	               status = h$value()[['status']],
	               status.message = h$value()[['statusMessage']])
	# check http status
	httpstatus <- as.numeric(output$status)
	if (httpstatus != 200) {
		print(output$status.message); break
		}
	else response <- fromJSON(output$data)
	return(response$result)
}

#' Returns a character vector corresponding to the currently available authorities.
#
#'
#' @return A list containing the currently available authorities.
#'
#' @references http://llama.mshri.on.ca/synergizer/translate/
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library('SynergizerR')
#' available_authorities()
#' }
#'
available_authorities <- function(){
	Sys.sleep(3) # In order to not have user host machine banned from using the service.
	header.field = c('Content-Type' = "application/json")
	curl <- getCurlHandle() 
	curlSetOpt(.opts = list(httpheader = header.field, verbose = FALSE), curl = curl) 
	t <- basicTextGatherer()
	h <- basicHeaderGatherer()

	body = curlPerform(url = "http://llama.mshri.on.ca/cgi/synergizer/serv",
	                   curl = curl,
	                   postfields = '{"method":"available_authorities","params":[],"id":0}',
	                   writefunction = t$update,
	                   headerfunction = h$update)	

	output <- list(data = t$value(),
	               status = h$value()[['status']],
	               status.message = h$value()[['statusMessage']])
	# check http status
	httpstatus <- as.numeric(output$status)
	if (httpstatus != 200) {
		print(output$status.message); break
		}
	else response <- fromJSON(output$data)
	return(response$result)
}



#' This method takes as parameter a single string, representing an authority, 
#' and returns a character vector corresponding to the currently available 
#' species for the chosen authority.
#'
#' @param authority A character containing any authoritative sources of identifier-mapping information.
#'
#' @return A vector containing  the currently available species for the chosen authority.
#'
#' @references http://llama.mshri.on.ca/synergizer/translate/
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library('SynergizerR')
#' available_species('ensembl')
#' }
#'
available_species <- function( authority = "ensembl" ){
	Sys.sleep(3) # In order to not have user host machine banned from using the service.
	header.field = c('Content-Type' = "application/json")
	curl <- getCurlHandle() 
	curlSetOpt(.opts = list(httpheader = header.field, verbose = FALSE), curl = curl) 
	t <- basicTextGatherer()
	h <- basicHeaderGatherer()

	body = curlPerform(url = "http://llama.mshri.on.ca/cgi/synergizer/serv",
	                   curl = curl,
	                   postfields = paste('{"method":"available_species","params":[', as.character(authority),
	                   '],"id":0}', sep="\""),
	                   writefunction = t$update,
	                   headerfunction = h$update)	
	
	output <- list(data = t$value(),
	               status = h$value()[['status']],
	               status.message = h$value()[['statusMessage']])
	# check http status
	httpstatus <- as.numeric(output$status)
	if (httpstatus != 200) {
		print(output$status.message); break
		}
	else response <- fromJSON(output$data)
	return(response$result)
}

#'Takes as parameters two strings, representing an authority and a species, and returns a character vector
#' corresponding to the currently available domain namespaces for the chosen authority and species.
#'
#' @param authority A character containing any authoritative sources of identifier-mapping information
#' @param species A character containing the Species. Note that the range of species supported 
#' depends on the choice of authority. Examples: Homo sapiens, Mus musculus.
#'
#' @return A vector containing the currently available domain namespaces for the chosen authority and species.
#'
#' @references http://llama.mshri.on.ca/synergizer/translate/
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library('SynergizerR')
#' available_domains('ensembl','homo sapiens')
#' }
#'
available_domains <- function( authority = "ensembl", species = "Homo sapiens") {
	Sys.sleep(3) # In order to not have user host machine banned from using the service.
	header.field = c('Content-Type' = "application/json")
	curl <- getCurlHandle() 
	curlSetOpt(.opts = list(httpheader = header.field, verbose = FALSE), curl = curl) 
	t <- basicTextGatherer()
	h <- basicHeaderGatherer()

	body = curlPerform(url = "http://llama.mshri.on.ca/cgi/synergizer/serv",
	                   curl = curl,
	                   postfields = paste('{"method":"available_domains","params":[', 
	                   as.character(authority), ",", as.character(species),
	                   '],"id":0}', sep="\""),
	                   writefunction = t$update,
	                   headerfunction = h$update)	
	
	output <- list(data = t$value(),
	               status = h$value()[['status']],
	               status.message = h$value()[['statusMessage']])
	# check http status
	httpstatus <- as.numeric(output$status)
	if (httpstatus != 200) {
		print(output$status.message); break
		}
	else response <- fromJSON(output$data)
	return(response$result)
}

#'Takes as parameters three strings, representing an authority, a species, and a domain namespace, 
#' and returns a character vector corresponding to the currently available range namespaces 
#' for the chosen authority, species, and domain namespace.
#'
#' @param authority A character containing any authoritative sources of identifier-mapping information.
#' @param species A character containing the Species. Note that the range of species supported 
#' depends on the choice of authority. Examples: Homo sapiens, Mus musculus.
#' @param domain This is the "namespace" (naming scheme) of the database identifiers the user wishes to translate.
#' Examples: embl, ipi
#'
#' @return A vector containing the currently available range namespaces 
#' for the chosen authority, species, and domain namespace.
#'
#' @references http://llama.mshri.on.ca/synergizer/translate/
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library('SynergizerR')
#' available_ranges('ensembl','homo sapiens','hgnc_symbol')
#' }
#' 
available_ranges <- function( authority = "ensembl", species = "Homo sapiens", domain ="hgnc_symbol") {
	Sys.sleep(3) # In order to not have user host machine banned from using the service.
	header.field = c('Content-Type' = "application/json")
	curl <- getCurlHandle() 
	curlSetOpt(.opts = list(httpheader = header.field, verbose = FALSE), curl = curl) 
	t <- basicTextGatherer()
	h <- basicHeaderGatherer()

	body = curlPerform(url = "http://llama.mshri.on.ca/cgi/synergizer/serv",
	                   curl = curl,
	                   postfields = paste('{"method":"available_ranges","params":[', 
	                   as.character(authority), ",", as.character(species), ",", as.character(domain),
	                   '],"id":0}', sep="\""),
	                   writefunction = t$update,
	                   headerfunction = h$update)	


	output <- list(data = t$value(),
	               status = h$value()[['status']],
	               status.message = h$value()[['statusMessage']])
	# check http status
	httpstatus <- as.numeric(output$status)
	if (httpstatus != 200) {
		print(output$status.message); break
		}
	else response <- fromJSON(output$data)
	return(response$result)
}
